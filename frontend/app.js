(function(){
  const cfg = window.APP_CONFIG;

  const el = {
    geneInput: document.getElementById('gene-input'),
    geneSuggestions: document.getElementById('gene-suggestions'),
    sourceSelect: document.getElementById('source-select'),
    patientInput: document.getElementById('patient-input'),
    plotBtn: document.getElementById('plot-btn'),
    plot: document.getElementById('plot'),
    status: document.getElementById('status'),
    geneInfo: document.querySelector('.gene-info'),
    geneName: document.querySelector('.gene-name'),
    geneDescription: document.querySelector('.gene-description')
  };

  // Populate sources dropdown
  cfg.sources.forEach((src, i) => {
    const opt = document.createElement('option');
    opt.value = src; opt.textContent = src.split('/').pop();
    if (i === 0) opt.selected = true;
    el.sourceSelect.appendChild(opt);
  });

  el.patientInput.value = cfg.defaultPatientId;

  let samplesMeta = null; // sample -> condition
  let currentData = null; // parsed CSV rows
  let geneToRows = new Map(); // gene_name -> rows
  let sampleColumns = []; // sample columns detected in file

  function setStatus(msg){ el.status.textContent = msg || ''; }

  function parseCsv(url){
    return new Promise((resolve, reject) => {
      Papa.parse(url, {
        header: true,
        dynamicTyping: true,
        download: true,
        skipEmptyLines: true,
        complete: res => resolve(res),
        error: err => reject(err)
      });
    });
  }

  async function loadSamples(){
    const res = await parseCsv(cfg.samplesCsv);
    const map = new Map();
    for (const row of res.data){
      if (!row.sample || !row.condition) continue;
      map.set(String(row.sample), String(row.condition));
    }
    return map;
  }

  function detectSampleColumns(rows){
    if (!rows || !rows.length) return [];
    const first = rows[0];
    const nonSample = new Set([cfg.columns.geneId, cfg.columns.geneName, cfg.columns.log2fc, cfg.columns.pvalue, cfg.columns.padj, "baseMean", "lfcSE", ""]);
    return Object.keys(first).filter(k => !nonSample.has(k));
  }

  function indexGenes(rows){
    const m = new Map();
    for (const r of rows){
      const g = (r[cfg.columns.geneName] || r[cfg.columns.geneId] || '').toString();
      if (!g) continue;
      if (!m.has(g)) m.set(g, []);
      m.get(g).push(r);
    }
    return m;
  }

  let allGenes = []; // Store all genes for filtering

  function updateSuggestions(filter = ''){
    el.geneSuggestions.innerHTML = '';
    let genes = allGenes;
    
    // Filter genes based on input
    if (filter) {
      genes = allGenes.filter(g => 
        g.toLowerCase().includes(filter.toLowerCase())
      );
    }
    
    // Limit to 50 options for performance
    genes = genes.slice(0, 50);
    
    for (const g of genes){
      const opt = document.createElement('option');
      opt.value = g; 
      el.geneSuggestions.appendChild(opt);
    }
  }

  function initializeGeneList(){
    // Sort genes alphabetically and filter out empty ones
    allGenes = Array.from(geneToRows.keys())
      .filter(g => g && g.trim())
      .sort((a, b) => a.localeCompare(b));
    updateSuggestions();
  }

  async function loadCurrentSource(){
    setStatus('Loading data...');
    const src = el.sourceSelect.value;
    const res = await parseCsv(src);
    currentData = res.data;
    sampleColumns = detectSampleColumns(currentData);
    geneToRows = indexGenes(currentData);
    initializeGeneList();
    setStatus(`Loaded ${currentData.length} rows from ${src.split('/').pop()}`);
  }

  function buildValuesForGene(gene){
    const rows = geneToRows.get(gene) || [];
    // If multiple rows map to same gene (duplicates), pick the row with max |log2FC| if available, else the first.
    let row = rows[0] || null;
    if (rows.length > 1){
      row = rows.slice().sort((a,b)=>Math.abs(+b[cfg.columns.log2fc]||0)-Math.abs(+a[cfg.columns.log2fc]||0))[0];
    }
    if (!row) return { controls: [], cases: [], myPoint: null };

    const controls = [];
    const cases = [];
    let myPt = null;
    const myPatientId = el.patientInput.value.trim();

    for (const s of sampleColumns){
      const sampleId = String(s);
      const cond = samplesMeta.get(sampleId);
      const y = Number(row[s]);
      
      if (!cond) continue; // skip unknown samples
      if (!Number.isFinite(y)) continue;
      
      const isMine = (sampleId === myPatientId);
      const pt = { x: cond, y, sample: sampleId, mine: isMine };
      
      if (cond === 'HC') {
        controls.push(pt);
      } else if (cond === 'D') {
        cases.push(pt);
        if (isMine) myPt = pt;
      }
    }
    
    return { controls, cases, myPoint: myPt };
  }

  function renderPlot(gene){
    const { controls, cases, myPoint } = buildValuesForGene(gene);
    if (!controls.length && !cases.length){
      setStatus('No data found for this gene in the selected file.');
      Plotly.purge(el.plot); return;
    }

    const otherCases = cases.filter(p=>!p.mine);
    const mine = myPoint ? [myPoint] : [];

    // Prism-style deterministic dot distribution
    const addPrismSpacing = (baseX, points) => {
      // Group points by Y value (rounded to avoid floating point issues)
      const yGroups = new Map();
      points.forEach((point, idx) => {
        const y = Math.round(point.y * 100) / 100; // Round to 2 decimals for grouping
        if (!yGroups.has(y)) yGroups.set(y, []);
        yGroups.get(y).push({point, originalIndex: idx, sampleId: point.sample});
      });
      
      const positions = new Array(points.length);
      const basePos = baseX === 'HC' ? 0 : 1;
      const prismSpread = 0.3; // Prism standard spread width
      
      yGroups.forEach((group) => {
        // Sort by sample ID for consistent ordering
        group.sort((a, b) => a.sampleId.localeCompare(b.sampleId));
        
        if (group.length === 1) {
          // Single point: center it
          positions[group[0].originalIndex] = basePos;
        } else {
          // Multiple points: evenly distribute like Prism
          group.forEach((item, idx) => {
            if (group.length === 2) {
              // Two points: spread symmetrically
              const offset = (idx === 0 ? -1 : 1) * prismSpread / 4;
              positions[item.originalIndex] = basePos + offset;
            } else {
              // Three or more points: even distribution
              const step = prismSpread / (group.length - 1);
              const offset = (idx * step) - (prismSpread / 2);
              positions[item.originalIndex] = basePos + offset;
            }
          });
        }
      });
      
      return positions;
    };

    const traceControls = {
      x: addPrismSpacing('HC', controls),
      y: controls.map(p=>p.y),
      mode: 'markers',
      type: 'scatter',
      marker: { 
        size: 10, 
        color: '#4472C4', 
        opacity: 0.8,
        line: { color: '#2E5AA0', width: 1 }
      },
      hoverinfo: 'skip', // do not show other patients' values
      name: 'Healthy Controls'
    };

    const traceOtherCases = {
      x: addPrismSpacing('D', otherCases),
      y: otherCases.map(p=>p.y),
      mode: 'markers',
      type: 'scatter',
      marker: { 
        size: 10, 
        color: '#4472C4', 
        opacity: 0.8,
        line: { color: '#2E5AA0', width: 1 }
      },
      hoverinfo: 'skip', // do not show other patients' values
      name: 'Other Cases'
    };

    const traceMine = {
      x: mine.map(()=> 1), // Position at D (1) centered for "My Patient"
      y: mine.map(p=>p.y),
      mode: 'markers',
      type: 'scatter',
      marker: { 
        size: 12, 
        color: '#8E44AD', 
        opacity: 1.0,
        line: { color: '#6A2C8A', width: 2 }
      },
      hovertemplate: '<b>Me</b><br>Sample: %{customdata[0]}<br>Value: %{y}<extra></extra>',
      customdata: mine.map(p=>[p.sample]),
      name: 'Me'
    };

    // Calculate medians for horizontal lines
    const getMedian = (values) => {
      if (values.length === 0) return null;
      const sorted = values.slice().sort((a,b) => a-b);
      const mid = Math.floor(sorted.length / 2);
      return sorted.length % 2 === 0 
        ? (sorted[mid - 1] + sorted[mid]) / 2 
        : sorted[mid];
    };
    
    const controlMedian = controls.length > 0 ? getMedian(controls.map(p => p.y)) : null;
    const caseMedian = cases.length > 0 ? getMedian(cases.map(p => p.y)) : null;
    


    // Add median lines
    const traces = [traceControls, traceOtherCases, traceMine].filter(t => t.x.length > 0);
    
    // Add horizontal median line for controls
    if (controlMedian !== null) {
      traces.push({
        x: [-0.2, 0.2],
        y: [controlMedian, controlMedian],
        mode: 'lines',
        type: 'scatter',
        line: { color: '#2C3E50', width: 3 },
        hoverinfo: 'none',
        showlegend: false,
        connectgaps: false
      });
    }
    
    // Add horizontal median line for cases  
    if (caseMedian !== null) {
      traces.push({
        x: [0.8, 1.2],
        y: [caseMedian, caseMedian],
        mode: 'lines',
        type: 'scatter',
        line: { color: '#2C3E50', width: 3 },
        hoverinfo: 'none',
        showlegend: false,
        connectgaps: false
      });
    }

    const layout = {
      paper_bgcolor: 'rgba(0,0,0,0)',
      plot_bgcolor: 'rgba(0,0,0,0)', 
      font: { color: '#ffffff', size: 12, family: 'Arial' },
      title: {
        text: `${gene}`,
        font: { size: 18, color: '#2C3E50', family: 'Arial Black' },
        x: 0.5,
        y: 0.95
      },
      margin: { l: 80, r: 40, t: 100, b: 80 },
      xaxis: {
        title: { text: '', font: { size: 16, color: '#2C3E50' } },
        tickvals: [0, 1],
        ticktext: ['HC', 'LC/ME'],
        tickfont: { size: 16, color: '#2C3E50', family: 'Arial Black' },
        gridcolor: 'rgba(0,0,0,0)', 
        zerolinecolor: 'rgba(0,0,0,0)',
        linecolor: '#2C3E50',
        linewidth: 2,
        ticklen: 5,
        showgrid: false,
        range: [-0.4, 1.4] // Add padding to show jittered points
      },
      yaxis: { 
        title: { 
          text: 'DESeq2 normalized counts', 
          font: { size: 16, color: '#2C3E50', family: 'Arial Black' },
          standoff: 20
        },
        tickfont: { size: 12, color: '#2C3E50' },
        gridcolor: 'rgba(0,0,0,0)', 
        zerolinecolor: 'rgba(0,0,0,0)',
        linecolor: '#2C3E50',
        linewidth: 2,
        showgrid: false,
        ticklen: 5
      },
      showlegend: false
    };

    Plotly.newPlot(el.plot, traces, layout, { displayModeBar: true, responsive: true });
  }

  // Typing effect for gene description
  function typeText(element, text, speed = 30) {
    return new Promise((resolve) => {
      element.innerHTML = '';
      let i = 0;
      
      function type() {
        if (i < text.length) {
          element.innerHTML += text.charAt(i);
          i++;
          setTimeout(type, speed);
        } else {
          // Add blinking cursor at the end
          element.innerHTML += '<span class="typing-cursor"></span>';
          resolve();
        }
      }
      
      type();
    });
  }

  // Generate placeholder description for genes
  function generateGeneDescription(geneName) {
    const descriptions = [
      `${geneName} is a protein-coding gene that plays a crucial role in cellular signaling pathways. Recent studies have shown its involvement in immune regulation and metabolic processes. This gene has been associated with various biological functions including cell adhesion, signal transduction, and transcriptional regulation.`,
      `The ${geneName} gene encodes a protein that is essential for normal cellular function. Research indicates this gene may be involved in development, immune response, or metabolic regulation. Expression levels of ${geneName} have been linked to several physiological processes and may serve as a biomarker in certain conditions.`,
      `${geneName} represents an important molecular target with significant biological implications. Studies suggest this gene plays a role in cellular homeostasis and may be involved in disease pathways. Further investigation of ${geneName} expression patterns could provide insights into therapeutic applications.`
    ];
    
    // Use gene name to deterministically select description
    const index = geneName.split('').reduce((acc, char) => acc + char.charCodeAt(0), 0) % descriptions.length;
    return descriptions[index];
  }

  async function showGeneInfo(geneName) {
    el.geneName.textContent = geneName;
    el.geneInfo.style.display = 'block';
    
    const description = generateGeneDescription(geneName);
    await typeText(el.geneDescription, description, 25);
  }

  async function init(){
    samplesMeta = await loadSamples();
    await loadCurrentSource();
    
    // Attach events
    el.sourceSelect.addEventListener('change', loadCurrentSource);
    
    // Real-time gene filtering as user types
    el.geneInput.addEventListener('input', ()=>{
      const filter = el.geneInput.value.trim();
      updateSuggestions(filter);
    });
    
    el.plotBtn.addEventListener('click', ()=>{
      const gene = el.geneInput.value.trim();
      if (!gene){ setStatus('Type a gene to plot.'); return; }
      
      // Validate patient ID is a case, not control
      const patientId = el.patientInput.value.trim();
      if (patientId && samplesMeta.has(patientId)) {
        const condition = samplesMeta.get(patientId);
        if (condition === 'HC') {
          setStatus('Error: Cannot select a healthy control patient ID. Please select a disease patient.');
          return;
        }
      }
      
      renderPlot(gene);
      showGeneInfo(gene);
    });
  }

  init().catch(err=>{
    console.error(err); setStatus('Failed to initialize. See console.');
  });
})();


