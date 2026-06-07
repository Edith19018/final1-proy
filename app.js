// =====================================================
// SimulaBolivia — app.js  (versión corregida)
// Correcciones: LU pivoteo parcial, tabla dif. divididas,
//               orden de convergencia, bug recalculo EDO
// =====================================================

window.addEventListener('load', () => {
  initHeroChart();
  initModule1Matrix();
  initModule3Data();
  document.getElementById('hamburger').addEventListener('click', () => {
    document.querySelector('.nav-links').classList.toggle('open');
  });
  document.getElementById('m2-func').addEventListener('change', () => {
    document.getElementById('m2-custom-group').style.display =
      document.getElementById('m2-func').value === 'custom' ? 'block' : 'none';
  });
});

// ======== HERO CHART ========
function initHeroChart() {
  const ctx = document.getElementById('heroChart').getContext('2d');
  const days = Array.from({length: 30}, (_, i) => i + 1);
  const prices = days.map(d => 8 + 0.45*d + 0.01*d*d + Math.sin(d*0.4)*0.5);
  new Chart(ctx, {
    type: 'line',
    data: {
      labels: days,
      datasets: [{
        label: 'Precio canasta (Bs)',
        data: prices,
        borderColor: '#f0b429',
        backgroundColor: 'rgba(240,180,41,0.08)',
        borderWidth: 2, fill: true, tension: 0.4, pointRadius: 0
      }]
    },
    options: {
      responsive: true,
      plugins: { legend: { labels: { color: '#7a8099', font: {family:'IBM Plex Mono', size:11} } } },
      scales: {
        x: { ticks: { color: '#4a5068', font:{size:10} }, grid: { color: 'rgba(255,255,255,0.04)' } },
        y: { ticks: { color: '#4a5068', font:{size:10} }, grid: { color: 'rgba(255,255,255,0.04)' } }
      }
    }
  });
}

// =====================================================
// MÓDULO 1: SISTEMAS DE ECUACIONES LINEALES
// =====================================================
const PRESETS = {
  normal:  { A: [[4,1,1],[2,5,2],[1,1,3]], b: [600,800,400],  label:'Distribución Normal' },
  bloqueo: { A: [[4,1,0],[2,5,2],[0,1,3]], b: [600,800,400],  label:'Ruta bloqueada' },
  crisis:  { A: [[4,1,1],[2,5,2],[1,1,3]], b: [900,1200,700], label:'Demanda Alta (crisis)' }
};

function initModule1Matrix() { loadPreset('normal'); }

function loadPreset(name) {
  const p = PRESETS[name];
  const matDiv = document.getElementById('m1-matrix');
  const vecDiv = document.getElementById('m1-vec-b');
  matDiv.innerHTML = '';
  vecDiv.innerHTML = '';
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      const inp = document.createElement('input');
      inp.type = 'number'; inp.id = `a${i}${j}`;
      inp.value = p.A[i][j]; inp.step = '0.1';
      matDiv.appendChild(inp);
    }
  }
  for (let i = 0; i < 3; i++) {
    const inp = document.createElement('input');
    inp.type = 'number'; inp.id = `b${i}`;
    inp.value = p.b[i]; inp.step = '0.1';
    vecDiv.appendChild(inp);
  }
}

function getMatrix() {
  const A = [], b = [];
  for (let i = 0; i < 3; i++) {
    A.push([]);
    for (let j = 0; j < 3; j++)
      A[i].push(parseFloat(document.getElementById(`a${i}${j}`).value) || 0);
    b.push(parseFloat(document.getElementById(`b${i}`).value) || 0);
  }
  return { A, b };
}

// --- CORRECCIÓN 1: LU con pivoteo parcial ---
function solveLU(Aorig, borig) {
  const n = Aorig.length;
  // Copia con índice de permutación
  const A = Aorig.map(r => [...r]);
  const b = [...borig];
  const perm = Array.from({length:n}, (_,i) => i); // vector de permutación
  const L = Array.from({length:n}, (_,i) => Array(n).fill(0));
  const steps = [];

  for (let k = 0; k < n; k++) {
    // Pivoteo parcial: buscar la fila con mayor |A[i][k]|
    let maxVal = Math.abs(A[k][k]), maxRow = k;
    for (let i = k+1; i < n; i++) {
      if (Math.abs(A[i][k]) > maxVal) { maxVal = Math.abs(A[i][k]); maxRow = i; }
    }
    if (maxRow !== k) {
      // Intercambiar filas en A, b y perm
      [A[k], A[maxRow]] = [A[maxRow], A[k]];
      [b[k], b[maxRow]] = [b[maxRow], b[k]];
      [perm[k], perm[maxRow]] = [perm[maxRow], perm[k]];
      // Intercambiar filas ya procesadas de L
      for (let j = 0; j < k; j++) [L[k][j], L[maxRow][j]] = [L[maxRow][j], L[k][j]];
      steps.push({iter: k+1, accion: `Pivote: fila ${k+1} ↔ fila ${maxRow+1}`, factor: '—'});
    }
    if (Math.abs(A[k][k]) < 1e-14) return null; // singular
    L[k][k] = 1;
    for (let i = k+1; i < n; i++) {
      const factor = A[i][k] / A[k][k];
      L[i][k] = factor;
      steps.push({iter: k+1, accion: `Eliminar fila ${i+1}`, factor: factor.toFixed(6)});
      for (let j = k; j < n; j++) A[i][j] -= factor * A[k][j];
    }
  }
  const U = A;
  // Sustitución hacia adelante Ly = b
  const y = Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    y[i] = b[i];
    for (let j = 0; j < i; j++) y[i] -= L[i][j] * y[j];
  }
  // Sustitución hacia atrás Ux = y
  const x = Array(n).fill(0);
  for (let i = n-1; i >= 0; i--) {
    x[i] = y[i];
    for (let j = i+1; j < n; j++) x[i] -= U[i][j] * x[j];
    x[i] /= U[i][i];
  }
  return { x, iters: steps, method: 'LU (con pivoteo parcial)', L, U };
}

// --- Jacobi ---
function solveJacobi(A, b, tol, maxIter) {
  const n = A.length;
  let x = Array(n).fill(0);
  const iters = [];
  for (let k = 0; k < maxIter; k++) {
    const xNew = Array(n).fill(0);
    for (let i = 0; i < n; i++) {
      let s = b[i];
      for (let j = 0; j < n; j++) if (j !== i) s -= A[i][j] * x[j];
      xNew[i] = s / A[i][i];
    }
    const err = Math.max(...xNew.map((v, i) => Math.abs(v - x[i])));
    iters.push({iter: k+1, x: xNew.map(v=>v.toFixed(4)), err: err.toExponential(3)});
    x = xNew;
    if (err < tol) break;
  }
  return { x, iters, method: 'Jacobi' };
}

// --- Gauss-Seidel ---
function solveGaussSeidel(A, b, tol, maxIter) {
  const n = A.length;
  let x = Array(n).fill(0);
  const iters = [];
  for (let k = 0; k < maxIter; k++) {
    const xOld = [...x];
    for (let i = 0; i < n; i++) {
      let s = b[i];
      for (let j = 0; j < n; j++) if (j !== i) s -= A[i][j] * x[j];
      x[i] = s / A[i][i];
    }
    const err = Math.max(...x.map((v, i) => Math.abs(v - xOld[i])));
    iters.push({iter: k+1, x: x.map(v=>v.toFixed(4)), err: err.toExponential(3)});
    if (err < tol) break;
  }
  return { x, iters, method: 'Gauss-Seidel' };
}

// --- SOR ---
function solveSOR(A, b, omega, tol, maxIter) {
  const n = A.length;
  let x = Array(n).fill(0);
  const iters = [];
  for (let k = 0; k < maxIter; k++) {
    const xOld = [...x];
    for (let i = 0; i < n; i++) {
      let s = b[i];
      for (let j = 0; j < n; j++) if (j !== i) s -= A[i][j] * x[j];
      x[i] = (1 - omega) * x[i] + omega * (s / A[i][i]);
    }
    const err = Math.max(...x.map((v, i) => Math.abs(v - xOld[i])));
    iters.push({iter: k+1, x: x.map(v=>v.toFixed(4)), err: err.toExponential(3)});
    if (err < tol) break;
  }
  return { x, iters, method: `SOR (ω=${omega})` };
}

// --- Gradiente Conjugado ---
function solveGradConj(A, b, tol, maxIter) {
  const n = A.length;
  const dot = (u, v) => u.reduce((s, _, i) => s + u[i]*v[i], 0);
  const Ax  = v => A.map(row => dot(row, v));
  const sub = (u, v) => u.map((x, i) => x - v[i]);
  let x = Array(n).fill(0);
  let r = sub(b, Ax(x));
  let p = [...r];
  const iters = [];
  for (let k = 0; k < maxIter; k++) {
    const Ap    = Ax(p);
    const rr    = dot(r, r);
    const alpha = rr / dot(p, Ap);
    x = x.map((xi, i) => xi + alpha * p[i]);
    const rNew  = sub(r, Ap.map(v => alpha * v));
    const beta  = dot(rNew, rNew) / rr;
    p = rNew.map((ri, i) => ri + beta * p[i]);
    const err = Math.sqrt(dot(rNew, rNew));
    iters.push({iter: k+1, x: x.map(v=>v.toFixed(4)), err: err.toExponential(3)});
    r = rNew;
    if (err < tol) break;
  }
  return { x, iters, method: 'Gradiente Conjugado' };
}

function conditionNumber(A) {
  const norm = M => Math.sqrt(M.flat().reduce((s, v) => s + v*v, 0));
  function inv3(M) {
    const det = M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
              - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
              + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);
    if (Math.abs(det) < 1e-12) return null;
    return [
      [(M[1][1]*M[2][2]-M[1][2]*M[2][1])/det, -(M[0][1]*M[2][2]-M[0][2]*M[2][1])/det, (M[0][1]*M[1][2]-M[0][2]*M[1][1])/det],
      [-(M[1][0]*M[2][2]-M[1][2]*M[2][0])/det, (M[0][0]*M[2][2]-M[0][2]*M[2][0])/det, -(M[0][0]*M[1][2]-M[0][2]*M[1][0])/det],
      [(M[1][0]*M[2][1]-M[1][1]*M[2][0])/det, -(M[0][0]*M[2][1]-M[0][1]*M[2][0])/det, (M[0][0]*M[1][1]-M[0][1]*M[1][0])/det]
    ];
  }
  const invA = inv3(A);
  return invA ? norm(A) * norm(invA) : Infinity;
}

let m1Chart = null;
function runModule1() {
  const { A, b } = getMatrix();
  const method   = document.getElementById('m1-method').value;
  const omega    = parseFloat(document.getElementById('m1-omega').value);
  const tol      = parseFloat(document.getElementById('m1-tol').value);
  const maxIter  = parseInt(document.getElementById('m1-maxiter').value);
  let res;
  switch (method) {
    case 'lu':         res = solveLU(A, b); break;
    case 'jacobi':     res = solveJacobi(A, b, tol, maxIter); break;
    case 'gauss_seidel': res = solveGaussSeidel(A, b, tol, maxIter); break;
    case 'sor':        res = solveSOR(A, b, omega, tol, maxIter); break;
    case 'grad_conj':  res = solveGradConj(A, b, tol, maxIter); break;
  }
  if (!res) { alert('El sistema no tiene solución única o la matriz es singular.'); return; }
  const x = res.x;
  const cond = conditionNumber(A);
  const zonas = ['Zona Norte', 'Zona Centro', 'Zona Sur'];
  const total = x.reduce((s,v) => s + parseFloat(v), 0);

  document.getElementById('m1-output').innerHTML = `
    <strong>Método:</strong> ${res.method}<br>
    <strong>Iteraciones / Pasos:</strong> ${res.iters.length}<br><br>
    ${zonas.map((z,i) => `→ <span class="result-highlight">${z}: ${parseFloat(x[i]).toFixed(2)} unidades</span>`).join('<br>')}
    <br><strong>Total distribuido:</strong> ${total.toFixed(2)} unidades<br><br>
    <strong>Número de condición κ(A) ≈</strong> ${isFinite(cond) ? cond.toFixed(2) : '∞'}
    ${cond > 100 ? '<br>⚠️ Sistema <strong>mal condicionado</strong> — sensible a pequeñas perturbaciones (rumores, bloqueos).'
                 : '<br>✓ Sistema bien condicionado — estable ante cambios pequeños en la demanda.'}
  `;

  // Tabla: LU muestra matrices L y U; iterativos muestran convergencia
  const tbl = document.getElementById('m1-table');
  if (method === 'lu' && res.L && res.U) {
    const fmt = v => v.toFixed(4);
    tbl.innerHTML = `
      <tr><th colspan="4">Matriz L (triangular inferior)</th></tr>
      <tr><th></th><th>col 1</th><th>col 2</th><th>col 3</th></tr>
      ${res.L.map((row,i)=>`<tr><th>fila ${i+1}</th>${row.map(v=>`<td>${fmt(v)}</td>`).join('')}</tr>`).join('')}
      <tr><th colspan="4" style="padding-top:1rem">Matriz U (triangular superior)</th></tr>
      <tr><th></th><th>col 1</th><th>col 2</th><th>col 3</th></tr>
      ${res.U.map((row,i)=>`<tr><th>fila ${i+1}</th>${row.map(v=>`<td>${fmt(v)}</td>`).join('')}</tr>`).join('')}
      <tr><th colspan="4" style="padding-top:1rem">Pasos de eliminación</th></tr>
      <tr><th>Paso k</th><th>Acción</th><th>Factor multiplicador</th><th></th></tr>
      ${res.iters.map(r=>`<tr><td>${r.iter}</td><td>${r.accion}</td><td>${r.factor}</td><td></td></tr>`).join('')}
    `;
  } else {
    const show = res.iters.slice(0, 15);
    tbl.innerHTML = `<tr><th>Iter</th><th>x₁</th><th>x₂</th><th>x₃</th><th>Error</th></tr>` +
      show.map(r => `<tr><td>${r.iter}</td><td>${r.x[0]}</td><td>${r.x[1]}</td><td>${r.x[2]}</td><td>${r.err}</td></tr>`).join('');
  }

  // Gráfico
  const ctx = document.getElementById('m1-chart').getContext('2d');
  if (m1Chart) m1Chart.destroy();
  m1Chart = new Chart(ctx, {
    type: 'bar',
    data: {
      labels: zonas,
      datasets: [
        { label: 'Cantidad asignada x', data: x.map(v=>parseFloat(parseFloat(v).toFixed(2))), backgroundColor:'rgba(240,180,41,0.7)', borderColor:'#f0b429', borderWidth:1 },
        { label: 'Demanda requerida b', data: b, backgroundColor:'rgba(255,255,255,0.1)', borderColor:'rgba(255,255,255,0.3)', borderWidth:1 }
      ]
    },
    options: {
      responsive: true,
      plugins: { legend: { labels: { color:'#7a8099' } } },
      scales: {
        x: { ticks:{ color:'#7a8099' }, grid:{ color:'rgba(255,255,255,0.04)' } },
        y: { ticks:{ color:'#7a8099' }, grid:{ color:'rgba(255,255,255,0.04)' } }
      }
    }
  });

  // Interpretación dinámica según resultado real
  const zonaMax = zonas[x.indexOf(Math.max(...x.map(v=>parseFloat(v))))];
  const zonaMin = zonas[x.indexOf(Math.min(...x.map(v=>parseFloat(v))))];
  document.getElementById('m1-interpretation').innerHTML = `
    <strong>Interpretación:</strong>
    Con el método <em>${res.method}</em>, la zona que recibe mayor cantidad de suministro es <strong>${zonaMax}</strong>
    y la que menos recibe es <strong>${zonaMin}</strong> (${parseFloat(x[zonas.indexOf(zonaMin)]).toFixed(2)} unidades).
    ${cond > 100
      ? ` El número de condición κ ≈ ${cond.toFixed(1)} es <strong>elevado</strong>: si la demanda varía un 1% — como ocurre en Bolivia cuando circula un rumor de escasez — la solución puede cambiar hasta un ${(cond*0.01*100).toFixed(0)}%. En un contexto de bloqueos de rutas, este sistema sería inestable y requeriría reasignación urgente de rutas alternativas.`
      : ` El número de condición κ ≈ ${cond.toFixed(1)} indica un sistema <strong>estable</strong>: la distribución resiste variaciones pequeñas en la demanda, como las fluctuaciones diarias típicas de los mercados de abastecimiento bolivianos.`}
    <br/><br/>En términos prácticos para Bolivia: si la zona ${zonaMin} recibe solo ${parseFloat(x[zonas.indexOf(zonaMin)]).toFixed(0)} unidades, podría quedar por debajo de la demanda mínima de sus surtidores o mercados, generando filas y descontento social. YPFB debería priorizar esa zona o habilitar una ruta alternativa desde otra planta.
  `;
  document.getElementById('m1-results').style.display = 'block';
}

// =====================================================
// MÓDULO 2: RAÍCES DE ECUACIONES
// =====================================================
function getF2() {
  const sel = document.getElementById('m2-func').value;
  const p1  = parseFloat(document.getElementById('m2-p1').value);
  const p2  = parseFloat(document.getElementById('m2-p2').value);
  switch(sel) {
    case 'gasto':   return x => p1 - p2 * Math.exp(0.03 * x);
    case 'reserva': return x => p1 * x - p2 * x * x;
    case 'opinion': return x => p1 * x * (1 - x) - p2 * x;
    case 'custom': {
      const expr = document.getElementById('m2-custom-fn').value;
      try { return new Function('x', `return ${expr};`); }
      catch(e) { alert('Función inválida'); return null; }
    }
  }
}

function biseccion(f, a, b, tol) {
  const iters = [];
  let fa = f(a);
  for (let k = 0; k < 200; k++) {
    const c = (a + b) / 2, fc = f(c);
    const err = Math.abs(b - a) / 2;
    iters.push({ iter:k+1, a:a.toFixed(6), b:b.toFixed(6), c:c.toFixed(6), fc:fc.toExponential(4), err });
    if (err < tol || Math.abs(fc) < 1e-12) return { root:c, iters };
    if (fa * fc < 0) { b = c; } else { a = c; fa = fc; }
  }
  return { root:(a+b)/2, iters };
}

function newtonRaphson(f, x0, tol) {
  const iters = [];
  let x = x0;
  const h = 1e-7;
  for (let k = 0; k < 100; k++) {
    const fx  = f(x);
    const dfx = (f(x+h) - f(x-h)) / (2*h);
    if (Math.abs(dfx) < 1e-14) break;
    const xNew = x - fx / dfx;
    const err  = Math.abs(xNew - x);
    iters.push({ iter:k+1, x:x.toFixed(8), fx:fx.toExponential(4), dfx:dfx.toExponential(4), err });
    x = xNew;
    if (err < tol) break;
  }
  return { root:x, iters };
}

function secante(f, x0, x1, tol) {
  const iters = [];
  for (let k = 0; k < 100; k++) {
    const f0 = f(x0), f1 = f(x1);
    if (Math.abs(f1 - f0) < 1e-14) break;
    const x2  = x1 - f1 * (x1 - x0) / (f1 - f0);
    const err = Math.abs(x2 - x1);
    iters.push({ iter:k+1, x0:x0.toFixed(6), x1:x1.toFixed(6), x2:x2.toFixed(6), err });
    x0 = x1; x1 = x2;
    if (err < tol) break;
  }
  return { root:x1, iters };
}

// CORRECCIÓN 3: Calcular orden de convergencia p ≈ log(e_{n+1}/e_n) / log(e_n/e_{n-1})
function calcOrdenConvergencia(iters) {
  const errs = iters.map(r => typeof r.err === 'number' ? r.err : parseFloat(r.err)).filter(e => e > 0);
  if (errs.length < 3) return null;
  const orders = [];
  for (let i = 1; i < errs.length - 1; i++) {
    const num = Math.log(errs[i+1] / errs[i]);
    const den = Math.log(errs[i]   / errs[i-1]);
    if (Math.abs(den) > 1e-10 && isFinite(num/den)) orders.push(num/den);
  }
  if (orders.length === 0) return null;
  // Promedio de los últimas 3 estimaciones (más estables)
  const tail = orders.slice(-3);
  return tail.reduce((s,v) => s+v, 0) / tail.length;
}

let m2Chart = null;
function runModule2() {
  const f      = getF2();
  if (!f) return;
  const method = document.getElementById('m2-method').value;
  const a      = parseFloat(document.getElementById('m2-a').value);
  const b      = parseFloat(document.getElementById('m2-b').value);
  const tol    = parseFloat(document.getElementById('m2-tol').value);
  let res;
  try {
    switch (method) {
      case 'biseccion':
        if (f(a)*f(b) >= 0) { alert('f(a) y f(b) deben tener signos opuestos para Bisección.'); return; }
        res = biseccion(f, a, b, tol); break;
      case 'newton':  res = newtonRaphson(f, a, tol); break;
      case 'secante': res = secante(f, a, b, tol); break;
    }
  } catch(e) { alert('Error evaluando la función: ' + e.message); return; }

  const root  = res.root;
  const orden = calcOrdenConvergencia(res.iters);
  const ordenTeorico = method === 'biseccion' ? 1 : method === 'newton' ? 2 : 1.618;
  const ordenStr = orden !== null
    ? `${orden.toFixed(3)} <span style="color:#7a8099;font-size:0.9em">(teórico: ${ordenTeorico})</span>`
    : 'insuficientes iteraciones para estimar';

  document.getElementById('m2-output').innerHTML = `
    <strong>Método:</strong> ${method}<br>
    <strong>Iteraciones:</strong> ${res.iters.length}<br>
    <strong>Raíz encontrada:</strong> <span class="result-highlight">x* = ${root.toFixed(8)}</span><br>
    <strong>f(x*) =</strong> ${f(root).toExponential(6)}<br>
    <strong>Orden de convergencia estimado p ≈</strong> ${ordenStr}
  `;

  // Gráfico con raíz marcada
  const xs = Array.from({length:200}, (_,i) => a + (b-a)*i/199);
  const ys = xs.map(x => { try { return f(x); } catch(e) { return null; } });
  const ctx = document.getElementById('m2-chart').getContext('2d');
  if (m2Chart) m2Chart.destroy();
  m2Chart = new Chart(ctx, {
    type: 'line',
    data: {
      labels: xs.map(v=>v.toFixed(2)),
      datasets: [
        { label:'f(x)', data:ys, borderColor:'#e05252', backgroundColor:'rgba(224,82,82,0.05)', borderWidth:2, fill:true, tension:0.3, pointRadius:0 },
        { label:`Raíz x*=${root.toFixed(4)}`, data:xs.map(x=>Math.abs(x-root)<(b-a)/150?0:null), borderColor:'#f0b429', borderWidth:0, pointRadius:xs.map(x=>Math.abs(x-root)<(b-a)/150?10:0), pointBackgroundColor:'#f0b429', showLine:false },
        { label:'y = 0', data:xs.map(()=>0), borderColor:'rgba(255,255,255,0.15)', borderWidth:1, pointRadius:0 }
      ]
    },
    options: {
      responsive:true,
      plugins:{ legend:{ labels:{ color:'#7a8099' } } },
      scales:{
        x:{ ticks:{ color:'#7a8099', maxTicksLimit:8 }, grid:{ color:'rgba(255,255,255,0.04)' } },
        y:{ ticks:{ color:'#7a8099' }, grid:{ color:'rgba(255,255,255,0.04)' } }
      }
    }
  });

  // Tabla de iteraciones con columna de orden de convergencia local
  const tbl      = document.getElementById('m2-table');
  const show     = res.iters.slice(0, 20);
  const errsArr  = show.map(r => typeof r.err === 'number' ? r.err : parseFloat(r.err));
  const ordenLocal = (i) => {
    if (i < 2) return '—';
    const num = Math.log(errsArr[i]   / errsArr[i-1]);
    const den = Math.log(errsArr[i-1] / errsArr[i-2]);
    const p   = Math.abs(den) > 1e-10 ? num/den : null;
    return (p !== null && isFinite(p)) ? p.toFixed(3) : '—';
  };

  if (method === 'biseccion') {
    tbl.innerHTML = `<tr><th>Iter</th><th>a</th><th>b</th><th>c=(a+b)/2</th><th>f(c)</th><th>Error</th><th>Orden local</th></tr>` +
      show.map((r,i) => `<tr><td>${r.iter}</td><td>${r.a}</td><td>${r.b}</td><td>${r.c}</td><td>${r.fc}</td><td>${r.err.toExponential(4)}</td><td>${ordenLocal(i)}</td></tr>`).join('');
  } else if (method === 'newton') {
    tbl.innerHTML = `<tr><th>Iter</th><th>xₙ</th><th>f(xₙ)</th><th>f\'(xₙ)</th><th>Error</th><th>Orden local</th></tr>` +
      show.map((r,i) => `<tr><td>${r.iter}</td><td>${r.x}</td><td>${r.fx}</td><td>${r.dfx}</td><td>${r.err.toExponential(4)}</td><td>${ordenLocal(i)}</td></tr>`).join('');
  } else {
    tbl.innerHTML = `<tr><th>Iter</th><th>x₀</th><th>x₁</th><th>x₂</th><th>Error</th><th>Orden local</th></tr>` +
      show.map((r,i) => `<tr><td>${r.iter}</td><td>${r.x0}</td><td>${r.x1}</td><td>${r.x2}</td><td>${r.err.toExponential(4)}</td><td>${ordenLocal(i)}</td></tr>`).join('');
  }

  const funcSel = document.getElementById('m2-func').value;
  let interp = `La raíz encontrada <strong>x* ≈ ${root.toFixed(4)}</strong> representa el umbral crítico buscado. `;
  if (funcSel === 'gasto')
    interp += `En el <strong>día ${root.toFixed(1)}</strong> del mes, el costo acumulado de la canasta básica iguala exactamente el ingreso familiar declarado. A partir de ese día, la familia boliviana entra en déficit: cada compra adicional implica endeudarse o reducir consumo. En contextos de alta inflación como el vivido en 2023-2024, este umbral se alcanzó cada vez más temprano en el mes.`;
  else if (funcSel === 'reserva')
    interp += `El equilibrio crítico de la reserva de combustible se alcanza en <strong>x* = ${root.toFixed(3)}</strong>. Por debajo de ese valor de tasa de reposición, el consumo supera sistemáticamente la entrada y la reserva se vacía inevitablemente. Este fue el escenario de varias plantas YPFB durante los picos de escasez de 2024, cuando las importaciones de diésel cayeron por debajo del punto de equilibrio.`;
  else
    interp += `Este es el punto de equilibrio social: el nivel de descontento x* donde la dinámica de movilización se estabiliza. ${orden && orden > 1.5 ? `El método ${method} converge <strong>cuadráticamente</strong> (orden ≈ ${orden.toFixed(2)}): cada iteración duplica los dígitos correctos, lo que lo hace muy eficiente para encontrar este umbral con pocos pasos.` : `El método ${method} converge <strong>linealmente</strong> (orden ≈ ${orden ? orden.toFixed(2) : '~1'}), requiriendo más iteraciones pero siendo más robusto ante condiciones iniciales alejadas.`}`;

  interp += ` <br/><br/><em>Comparativa de convergencia: Bisección garantiza la raíz pero necesita más iteraciones; Newton-Raphson es el más rápido si la función tiene derivada bien definida cerca de la raíz; Secante es un buen compromiso cuando calcular la derivada analítica es costoso.</em>`;

  document.getElementById('m2-interpretation').innerHTML = `<strong>Interpretación:</strong> ${interp}`;
  document.getElementById('m2-results').style.display = 'block';
}

// =====================================================
// MÓDULO 3: INTERPOLACIÓN
// =====================================================
const DEFAULT_DATA = [
  {x:1,y:8},{x:5,y:10},{x:10,y:13},{x:15,y:16},{x:20,y:19},{x:30,y:22}
];

function initModule3Data() { resetDataRows(); }

function resetDataRows() {
  const cont = document.getElementById('m3-data-rows');
  cont.innerHTML = '';
  DEFAULT_DATA.forEach(p => addDataRow(p.x, p.y));
}

function addDataRow(xv='', yv='') {
  const cont = document.getElementById('m3-data-rows');
  const row  = document.createElement('div');
  row.className = 'data-row';
  row.innerHTML = `
    <label>Día:</label><input type="number" class="x-val" value="${xv}" placeholder="Día"/>
    <label>Precio (Bs):</label><input type="number" class="y-val" value="${yv}" step="0.01" placeholder="Precio"/>
    <button class="btn-del" onclick="this.parentElement.remove()">✕</button>
  `;
  cont.appendChild(row);
}

function getDataPoints() {
  const rows = document.querySelectorAll('#m3-data-rows .data-row');
  const pts  = [];
  rows.forEach(r => {
    const x = parseFloat(r.querySelector('.x-val').value);
    const y = parseFloat(r.querySelector('.y-val').value);
    if (!isNaN(x) && !isNaN(y)) pts.push({x, y});
  });
  pts.sort((a,b) => a.x - b.x);
  return pts;
}

function lagrange(pts, xEval) {
  const n = pts.length;
  let result = 0;
  for (let i = 0; i < n; i++) {
    let L = 1;
    for (let j = 0; j < n; j++) if (j !== i) L *= (xEval - pts[j].x) / (pts[i].x - pts[j].x);
    result += pts[i].y * L;
  }
  return result;
}

// CORRECCIÓN 2: newtonDivDiff retorna la tabla completa para renderizar
function newtonDivDiff(pts, xEval) {
  const n  = pts.length;
  // Construir tabla de diferencias divididas
  const dd = pts.map(p => [p.y]);
  for (let j = 1; j < n; j++) {
    for (let i = 0; i < n - j; i++) {
      dd[i].push((dd[i+1][j-1] - dd[i][j-1]) / (pts[i+j].x - pts[i].x));
    }
  }
  const coefs = dd.map(row => row[0]); // diagonal principal = coeficientes
  // Evaluar polinomio de Newton
  let result = coefs[0];
  let term   = 1;
  for (let i = 1; i < n; i++) {
    term   *= (xEval - pts[i-1].x);
    result += coefs[i] * term;
  }
  return { result, coefs, dd, n };
}

function cubicSplines(pts, xEval) {
  const n = pts.length;
  const h = [];
  for (let i = 0; i < n-1; i++) h.push(pts[i+1].x - pts[i].x);
  const size = n - 2;
  if (size <= 0) return lagrange(pts, xEval);
  // Verificar que no haya h ≈ 0 (puntos muy cercanos)
  if (h.some(hi => Math.abs(hi) < 1e-10)) return lagrange(pts, xEval);
  const diag = [], lower = [], upper = [], rhs = [];
  for (let i = 0; i < size; i++) {
    diag.push(2 * (h[i] + h[i+1]));
    if (i > 0) lower.push(h[i]);
    if (i < size-1) upper.push(h[i+1]);
    rhs.push(6 * ((pts[i+2].y - pts[i+1].y)/h[i+1] - (pts[i+1].y - pts[i].y)/h[i]));
  }
  const m = Array(size).fill(0);
  const d = [...diag], r = [...rhs];
  for (let i = 1; i < size; i++) {
    const factor = lower[i-1] / d[i-1];
    d[i] -= factor * upper[i-1];
    r[i] -= factor * r[i-1];
  }
  m[size-1] = r[size-1] / d[size-1];
  for (let i = size-2; i >= 0; i--) m[i] = (r[i] - upper[i]*m[i+1]) / d[i];
  const M = [0, ...m, 0];
  let seg = n - 2;
  for (let i = 0; i < n-1; i++) { if (xEval <= pts[i+1].x) { seg = i; break; } }
  const xi = pts[seg].x, xi1 = pts[seg+1].x, hi = h[seg];
  return M[seg]/(6*hi)*Math.pow(xi1-xEval,3) + M[seg+1]/(6*hi)*Math.pow(xEval-xi,3)
       + (pts[seg].y/hi - M[seg]*hi/6)*(xi1-xEval) + (pts[seg+1].y/hi - M[seg+1]*hi/6)*(xEval-xi);
}

let m3Chart = null;
function runModule3() {
  const pts = getDataPoints();
  if (pts.length < 3) { alert('Se necesitan al menos 3 puntos.'); return; }
  const method = document.getElementById('m3-method').value;
  const xEval  = parseFloat(document.getElementById('m3-xeval').value);
  const xMin   = pts[0].x, xMax = pts[pts.length-1].x;

  // Para Newton guardamos el resultado completo
  let newtonResult = null;
  const evalFn = x => {
    switch(method) {
      case 'lagrange': return lagrange(pts, x);
      case 'newton':
        if (!newtonResult || newtonResult._xEval !== x)
          newtonResult = { ...newtonDivDiff(pts, x), _xEval: x };
        return newtonResult.result;
      case 'splines': return cubicSplines(pts, x);
    }
  };

  // Obtener resultado completo de Newton para xEval (para la tabla)
  const newtonFull = method === 'newton' ? newtonDivDiff(pts, xEval) : null;

  const xCurve = Array.from({length:120}, (_,i) => xMin + (xMax-xMin)*i/119);
  const yCurve = xCurve.map(x => { try { return evalFn(x); } catch(e) { return null; } });
  const yEstim = evalFn(xEval);

  document.getElementById('m3-output').innerHTML = `
    <strong>Método:</strong> ${method}<br>
    <strong>Día evaluado:</strong> ${xEval}<br>
    <strong>Precio estimado:</strong> <span class="result-highlight">Bs ${yEstim.toFixed(2)}</span><br>
    <strong>Rango de datos:</strong> Día ${xMin} → Día ${xMax}
    ${xEval < xMin || xEval > xMax ? '<br>⚠️ Fuera del rango — <strong>extrapolación</strong>: menor confiabilidad.' : ''}
  `;

  const ctx = document.getElementById('m3-chart').getContext('2d');
  if (m3Chart) m3Chart.destroy();
  m3Chart = new Chart(ctx, {
    type: 'line',
    data: {
      labels: xCurve.map(v=>v.toFixed(1)),
      datasets: [
        { label:`Curva (${method})`, data:yCurve, borderColor:'#4fc3f7', backgroundColor:'rgba(79,195,247,0.06)', borderWidth:2, fill:true, tension:0.3, pointRadius:0 },
        { label:'Datos originales', data:xCurve.map(xc => { const p=pts.find(p=>Math.abs(p.x-xc)<0.25); return p?p.y:null; }), borderColor:'#f0b429', borderWidth:0, pointRadius:6, pointBackgroundColor:'#f0b429', showLine:false },
        { label:`Día ${xEval} → Bs ${yEstim.toFixed(2)}`, data:xCurve.map(xc=>Math.abs(xc-xEval)<(xMax-xMin)/80?yEstim:null), borderColor:'#e05252', borderWidth:0, pointRadius:xCurve.map(xc=>Math.abs(xc-xEval)<(xMax-xMin)/80?9:0), pointBackgroundColor:'#e05252', showLine:false }
      ]
    },
    options: {
      responsive:true,
      plugins:{ legend:{ labels:{ color:'#7a8099' } } },
      scales:{
        x:{ ticks:{ color:'#7a8099', maxTicksLimit:8 }, grid:{ color:'rgba(255,255,255,0.04)' } },
        y:{ ticks:{ color:'#7a8099' }, grid:{ color:'rgba(255,255,255,0.04)' } }
      }
    }
  });

  // CORRECCIÓN 2: Mostrar tabla de diferencias divididas si es Newton
  const tbl = document.getElementById('m3-table');
  if (method === 'newton' && newtonFull) {
    const { dd, coefs, n } = newtonFull;
    // Cabecera dinámica
    let head = '<tr><th>xᵢ</th><th>f[xᵢ]</th>';
    for (let j = 1; j < n; j++) head += `<th>Orden ${j}</th>`;
    head += '</tr>';
    let body = '';
    for (let i = 0; i < n; i++) {
      body += `<tr><td>${pts[i].x}</td>`;
      for (let j = 0; j < n - i; j++) {
        body += `<td>${dd[i][j] !== undefined ? dd[i][j].toFixed(5) : ''}</td>`;
      }
      // Rellenar celdas vacías
      for (let j = n-i; j < n; j++) body += '<td></td>';
      body += '</tr>';
    }
    // Coeficientes
    body += `<tr><td colspan="${n+1}" style="padding-top:.75rem"><strong>Coeficientes del polinomio de Newton:</strong> ${coefs.map((c,i)=>`c${i}=${c.toFixed(5)}`).join(' | ')}</td></tr>`;
    tbl.innerHTML = `<caption style="caption-side:top;padding:.5rem;color:#7a8099;font-size:.82rem">Tabla de diferencias divididas</caption>${head}${body}`;
  } else {
    // Tabla de valores interpolados
    const step = Math.ceil((xMax-xMin)/9);
    const evalPts = Array.from({length:10}, (_,i) => xMin + i*step).filter(x=>x<=xMax);
    tbl.innerHTML = `<tr><th>Día</th><th>Precio estimado (Bs)</th><th>Dato original</th><th>Error vs dato</th></tr>` +
      evalPts.map(xc => {
        const orig = pts.find(p => p.x === xc);
        const est  = evalFn(xc);
        const err  = orig ? Math.abs(est - orig.y).toFixed(4) : '—';
        return `<tr><td>${xc}</td><td>${est.toFixed(2)}</td><td>${orig ? orig.y.toFixed(2) : '—'}</td><td>${err}</td></tr>`;
      }).join('');
  }

  const pts3 = getDataPoints();
  const precioMin = Math.min(...pts3.map(p=>p.y));
  const precioMax = Math.max(...pts3.map(p=>p.y));
  const incr = ((precioMax - precioMin) / precioMin * 100).toFixed(1);
  document.getElementById('m3-interpretation').innerHTML = `
    <strong>Interpretación:</strong> El precio estimado para el día <strong>${xEval}</strong> es <strong>Bs ${yEstim.toFixed(2)}</strong>/kg.
    <br/><br/>
    Durante el período analizado, el precio pasó de <strong>Bs ${precioMin.toFixed(2)}</strong> a <strong>Bs ${precioMax.toFixed(2)}</strong>, un incremento del <strong>${incr}%</strong>.
    En Bolivia, alzas de este tipo en la papa y otros básicos fueron frecuentes en mercados como la Rodríguez (La Paz) o el mercado 25 de Mayo (Cochabamba) durante los períodos de bloqueo de 2024.
    <br/><br/>
    ${method === 'splines'
      ? '<strong>Splines cúbicos:</strong> producen la curva más suave y realista al garantizar continuidad hasta la segunda derivada. No oscilan artificialmente entre datos. Recomendado para análisis de precios con tendencia gradual — refleja bien la realidad de mercados que suben sostenidamente sin bajas bruscas entre semana.'
      : method === 'lagrange'
      ? '<strong>Lagrange:</strong> polinomio de grado n−1 que pasa exactamente por todos los puntos. Con 6 datos produce un polinomio de grado 5, suficientemente suave. Puede presentar el <em>fenómeno de Runge</em> con datos muy dispersos o en los extremos del intervalo — lo que en términos económicos equivaldría a predecir incorrectamente que el precio bajó entre dos semanas consecutivas.'
      : '<strong>Newton (diferencias divididas):</strong> matemáticamente equivalente a Lagrange pero construido incrementalmente. Ventaja práctica: si se añade un nuevo registro de precios de feria, solo se agrega un término sin recalcular todo. La tabla triangular revela el <em>orden de variación</em>: diferencias de primer orden indican la tasa de cambio diaria del precio.'}
  `;
  document.getElementById('m3-results').style.display = 'block';
}

// =====================================================
// MÓDULO 4: INTEGRACIÓN NUMÉRICA
// =====================================================
function getCostFn4() {
  const sel = document.getElementById('m4-func').value;
  const C0  = parseFloat(document.getElementById('m4-c0').value);
  const r   = parseFloat(document.getElementById('m4-r').value);
  switch(sel) {
    case 'lineal':      return t => C0 + r * t;
    case 'exponencial': return t => C0 * Math.exp(r * 0.01 * t);
    case 'escalon': {
      const pts = getDataPoints();
      return t => pts.length >= 3
        ? cubicSplines(pts, Math.min(Math.max(t, pts[0].x), pts[pts.length-1].x))
        : C0 + r * t;
    }
  }
}

function trapecio(f, a, b, n) {
  const h = (b-a)/n;
  let s = f(a) + f(b);
  for (let i=1; i<n; i++) s += 2*f(a+i*h);
  return s * h / 2;
}

function simpson13(f, a, b, n) {
  if (n % 2 !== 0) n++;
  const h = (b-a)/n;
  let s = f(a) + f(b);
  for (let i=1; i<n; i++) s += (i%2===0?2:4)*f(a+i*h);
  return s * h / 3;
}

function simpson38(f, a, b, n) {
  if (n % 3 !== 0) n = Math.ceil(n/3)*3;
  const h = (b-a)/n;
  let s = f(a) + f(b);
  for (let i=1; i<n; i++) s += (i%3===0?2:3)*f(a+i*h);
  return s * 3*h / 8;
}

let m4Chart = null;
function runModule4() {
  const C0      = parseFloat(document.getElementById('m4-c0').value);
  const T       = parseFloat(document.getElementById('m4-T').value);
  const n       = parseInt(document.getElementById('m4-n').value);
  const ingreso = parseFloat(document.getElementById('m4-ingreso').value);
  const f = getCostFn4();
  const trap   = trapecio(f, 0, T, n);
  const simp13 = simpson13(f, 0, T, n);
  const simp38 = simpson38(f, 0, T, n);
  const gastoEstable = C0 * T;
  const perdida = trap - gastoEstable;
  const pctTrap = (perdida / ingreso * 100).toFixed(1);
  // Error relativo entre métodos (Simpson 1/3 como referencia)
  const errTrap   = Math.abs(trap - simp13) / Math.abs(simp13) * 100;
  const errSimp38 = Math.abs(simp38 - simp13) / Math.abs(simp13) * 100;

  document.getElementById('m4-output').innerHTML = `
    <strong>Gasto acumulado — Trapecio:</strong>    <span class="result-highlight">Bs ${trap.toFixed(2)}</span><br>
    <strong>Gasto acumulado — Simpson 1/3:</strong> <span class="result-highlight">Bs ${simp13.toFixed(2)}</span><br>
    <strong>Gasto acumulado — Simpson 3/8:</strong> <span class="result-highlight">Bs ${simp38.toFixed(2)}</span><br>
    <strong>Gasto sin inflación:</strong> Bs ${gastoEstable.toFixed(2)}<br>
    <strong>Pérdida del poder adquisitivo:</strong> <span class="result-highlight">Bs ${perdida.toFixed(2)} (${pctTrap}% del ingreso mensual)</span><br>
    ${parseFloat(pctTrap) > 20 ? '<br>⚠️ La familia destina más del 20% extra de su ingreso solo por inflación.' : ''}
  `;

  const ts    = Array.from({length:T+1}, (_,i)=>i);
  const costs = ts.map(t=>f(t));
  const cumTrap = ts.map(t => t===0 ? 0 : trapecio(f,0,t,Math.max(2,t)));
  const ctx = document.getElementById('m4-chart').getContext('2d');
  if (m4Chart) m4Chart.destroy();
  m4Chart = new Chart(ctx, {
    type:'line',
    data:{
      labels: ts,
      datasets:[
        { label:'Costo diario C(t) [Bs]', data:costs, borderColor:'#66bb6a', backgroundColor:'rgba(102,187,106,0.07)', borderWidth:2, fill:true, tension:0.3, pointRadius:0, yAxisID:'y1' },
        { label:'Gasto acumulado (Trapecio)', data:cumTrap, borderColor:'#f0b429', backgroundColor:'rgba(240,180,41,0.07)', borderWidth:2, fill:true, tension:0.3, pointRadius:0, yAxisID:'y2' },
        { label:'Sin inflación (ref.)', data:ts.map(t=>C0*t), borderColor:'rgba(255,255,255,0.2)', borderWidth:1, borderDash:[4,4], pointRadius:0, yAxisID:'y2' }
      ]
    },
    options:{
      responsive:true,
      plugins:{ legend:{ labels:{ color:'#7a8099' } } },
      scales:{
        x:{ ticks:{ color:'#7a8099' }, grid:{ color:'rgba(255,255,255,0.04)' } },
        y1:{ type:'linear', position:'left',  ticks:{ color:'#66bb6a' }, grid:{ color:'rgba(255,255,255,0.04)' } },
        y2:{ type:'linear', position:'right', ticks:{ color:'#f0b429' }, grid:{ drawOnChartArea:false } }
      }
    }
  });

  document.getElementById('m4-table').innerHTML = `
    <tr><th>Método</th><th>Gasto acumulado (Bs)</th><th>Pérdida vs sin inflación (Bs)</th><th>% del ingreso</th><th>Error relativo vs S1/3</th></tr>
    <tr><td>Trapecio</td><td>${trap.toFixed(4)}</td><td>${(trap-gastoEstable).toFixed(4)}</td><td>${((trap-gastoEstable)/ingreso*100).toFixed(3)}%</td><td>${errTrap.toFixed(4)}%</td></tr>
    <tr><td>Simpson 1/3</td><td>${simp13.toFixed(4)}</td><td>${(simp13-gastoEstable).toFixed(4)}</td><td>${((simp13-gastoEstable)/ingreso*100).toFixed(3)}%</td><td>— (referencia)</td></tr>
    <tr><td>Simpson 3/8</td><td>${simp38.toFixed(4)}</td><td>${(simp38-gastoEstable).toFixed(4)}</td><td>${((simp38-gastoEstable)/ingreso*100).toFixed(3)}%</td><td>${errSimp38.toFixed(4)}%</td></tr>
  `;

  // Interpretación dinámica
  const masPrec = errTrap < errSimp38 ? 'Trapecio' : 'Simpson 3/8';
  document.getElementById('m4-interpretation').innerHTML = `
    <strong>Interpretación:</strong>
    En <strong>${T} días</strong> con precios crecientes, la familia gastó <strong>Bs ${trap.toFixed(2)}</strong> en su canasta básica.
    Sin inflación habría gastado <strong>Bs ${gastoEstable.toFixed(2)}</strong> — la diferencia de <strong>Bs ${perdida.toFixed(2)}</strong>
    representa el <strong>${pctTrap}%</strong> del ingreso mensual perdido por causa de la inflación.
    <br/><br/>
    En el contexto boliviano, con un ingreso de referencia de Bs ${ingreso.toFixed(0)} y precios que suben desde el primer día del mes,
    ${perdida > ingreso * 0.2
      ? 'esta pérdida supera el <strong>20% del ingreso mensual</strong>: una situación crítica que obliga a las familias del área periurbana a reducir la cantidad o calidad de alimentos adquiridos. Este escenario es coherente con lo reportado en El Alto y Oruro durante los picos de escasez de 2024.'
      : 'la pérdida es significativa pero manejable si el ingreso familiar supera el salario mínimo nacional (Bs 2.362). Sin embargo, para las familias que dependen del comercio informal con ingresos cercanos al mínimo, incluso una pérdida del ' + pctTrap + '% puede implicar reducir la cantidad de alimentos.'}
    <br/><br/>
    <strong>Comparativa de métodos:</strong> Simpson 1/3 tiene error O(h⁴) y es más preciso que el Trapecio O(h²) para la misma cantidad de subintervalos.
    En esta simulación, ${masPrec === 'Simpson 3/8' ? 'Simpson 3/8' : 'el Trapecio'} se acerca más a Simpson 1/3 con diferencia relativa del ${Math.min(errTrap,errSimp38).toFixed(4)}%.
    Para n = ${n} subintervalos, todos los métodos convergen bien — la diferencia entre ellos es matemáticamente pequeña pero conceptualmente importante: Simpson 1/3 es el estándar recomendado.
  `;
  document.getElementById('m4-results').style.display = 'block';
}

// =====================================================
// MÓDULO 5: ECUACIONES DIFERENCIALES
// =====================================================
function loadEDOPreset(name) {
  const presets = {
    normal:  {r0:10000, entrada:1000, consumo:1000, panico:0.005, dias:30,  h:0.5},
    bloqueo: {r0:10000, entrada:0,    consumo:1000, panico:0.01,  dias:20,  h:0.5},
    panico:  {r0:10000, entrada:600,  consumo:800,  panico:0.05,  dias:25,  h:0.25}
  };
  const p = presets[name];
  document.getElementById('m5-r0').value      = p.r0;
  document.getElementById('m5-entrada').value = p.entrada;
  document.getElementById('m5-consumo').value = p.consumo;
  document.getElementById('m5-panico').value  = p.panico;
  document.getElementById('m5-dias').value    = p.dias;
  document.getElementById('m5-h').value       = p.h;
}

function solveEDO(method, R0, entrada, consumo, panico, dias, h) {
  const f     = (t, R) => entrada - consumo * (1 + panico * t);
  const steps = Math.ceil(dias / h);
  const ts    = [0], Rs = [R0];
  for (let i = 0; i < steps; i++) {
    const t = ts[ts.length-1], R = Rs[Rs.length-1];
    if (R <= 0) { ts.push(+(t+h).toFixed(4)); Rs.push(0); continue; }
    let Rnew;
    switch(method) {
      case 'euler':
        Rnew = R + h * f(t, R);
        break;
      case 'heun': {
        const k1 = f(t,   R);
        const k2 = f(t+h, R + h*k1);
        Rnew = R + h*(k1+k2)/2;
        break;
      }
      case 'rk4': {
        const k1 = f(t,     R);
        const k2 = f(t+h/2, R + h/2*k1);
        const k3 = f(t+h/2, R + h/2*k2);
        const k4 = f(t+h,   R + h*k3);
        Rnew = R + h*(k1 + 2*k2 + 2*k3 + k4)/6;
        break;
      }
    }
    ts.push(+(t+h).toFixed(4));
    Rs.push(Math.max(0, Rnew));
  }
  return { ts, Rs };
}

let m5Chart = null;
function runModule5() {
  const method  = document.getElementById('m5-method').value;
  const R0      = parseFloat(document.getElementById('m5-r0').value);
  const entrada = parseFloat(document.getElementById('m5-entrada').value);
  const consumo = parseFloat(document.getElementById('m5-consumo').value);
  const panico  = parseFloat(document.getElementById('m5-panico').value);
  const dias    = parseFloat(document.getElementById('m5-dias').value);
  const h       = parseFloat(document.getElementById('m5-h').value);
  const CRITICO = R0 * 0.15;

  const colors  = { euler:'#e05252', heun:'#4fc3f7', rk4:'#66bb6a' };
  const methods = method === 'todos' ? ['euler','heun','rk4'] : [method];

  // CORRECCIÓN 4: calcular TODOS los resultados UNA sola vez
  const resultados = {};
  methods.forEach(m => {
    resultados[m] = solveEDO(m, R0, entrada, consumo, panico, dias, h);
  });
  const mainKey    = methods[0];
  const { ts, Rs } = resultados[mainKey];

  // Gráfico
  const ctx = document.getElementById('m5-chart').getContext('2d');
  if (m5Chart) m5Chart.destroy();
  const datasets = methods.map(m => ({
    label: m.toUpperCase(),
    data: resultados[m].ts.map((t,i) => ({x:t, y:resultados[m].Rs[i]})),
    borderColor: colors[m],
    backgroundColor: colors[m]+'20',
    borderWidth: 2, fill: false, tension: 0.1, pointRadius: 0
  }));
  datasets.push({
    label: `Nivel crítico (${CRITICO.toFixed(0)} L)`,
    data: [{x:0,y:CRITICO},{x:dias,y:CRITICO}],
    borderColor:'rgba(240,180,41,0.6)',
    borderWidth:1, borderDash:[6,4], pointRadius:0
  });

  m5Chart = new Chart(ctx, {
    type:'line',
    data:{ datasets },
    options:{
      responsive:true,
      parsing:{ xAxisKey:'x', yAxisKey:'y' },
      plugins:{ legend:{ labels:{ color:'#7a8099' } } },
      scales:{
        x:{ type:'linear', ticks:{color:'#7a8099'}, grid:{color:'rgba(255,255,255,0.04)'}, title:{display:true,text:'Días',color:'#7a8099'} },
        y:{ ticks:{color:'#7a8099'}, grid:{color:'rgba(255,255,255,0.04)'}, title:{display:true,text:'Reserva (litros)',color:'#7a8099'} }
      }
    }
  });

  // Día crítico con el método principal
  const idxCrit  = Rs.findIndex(r => r <= CRITICO);
  const diaCrit  = idxCrit >= 0 ? ts[idxCrit].toFixed(1) : 'No alcanzado';
  const idxVacio = Rs.findIndex(r => r <= 0);
  const diaVacio = idxVacio >= 0 ? ts[idxVacio].toFixed(1) : 'No alcanzado';

  document.getElementById('m5-output').innerHTML = `
    <strong>Reserva inicial:</strong> ${R0.toLocaleString()} litros<br>
    <strong>Déficit diario base:</strong> ${Math.max(0, consumo - entrada).toFixed(0)} L/día<br>
    <strong>Día nivel crítico (15%):</strong> <span class="result-highlight">Día ${diaCrit}</span><br>
    <strong>Día reserva = 0:</strong>         <span class="result-highlight">Día ${diaVacio}</span><br>
    <strong>Reserva al final (día ${dias}):</strong> ${Rs[Rs.length-1].toFixed(0)} litros
  `;

  // CORRECCIÓN 4: tabla usa resultados ya calculados — sin recalcular
  const tbl  = document.getElementById('m5-table');
  const step = Math.max(1, Math.floor(ts.length / 15));
  tbl.innerHTML = `<tr><th>Día</th>${methods.map(m=>`<th>${m.toUpperCase()} R(t) [L]</th>`).join('')}${methods.length>1?'<th>Δ Euler−RK4</th>':''}</tr>`;
  for (let i = 0; i < ts.length; i += step) {
    const vals = methods.map(m => resultados[m].Rs[i] ?? 0);
    const diff = methods.length > 1 ? Math.abs(vals[0] - vals[vals.length-1]).toFixed(2) : null;
    tbl.innerHTML += `<tr><td>${ts[i].toFixed(2)}</td>${vals.map(v=>`<td>${v.toFixed(1)}</td>`).join('')}${diff!==null?`<td>${diff}</td>`:''}</tr>`;
  }

  // Interpretación dinámica
  const errEulerRK4 = methods.length > 1
    ? Math.abs(resultados['euler'].Rs[resultados['euler'].Rs.length-1] - resultados['rk4'].Rs[resultados['rk4'].Rs.length-1]).toFixed(1)
    : null;

  const escenarioTexto = entrada === 0
    ? `<strong>Escenario de bloqueo total:</strong> sin ningún reabastecimiento, la reserva cae a razón del consumo puro. Este fue el caso de algunas plantas YPFB cuando bloqueos de carretera cortaron completamente el suministro por tierra. En esas condiciones, incluso una reserva amplia puede agotarse en días.`
    : panico > 0.03
    ? `<strong>Escenario de pánico:</strong> el factor κ = ${panico} hace que el consumo crezca rápidamente con el tiempo, acelerando el vaciado de forma no lineal. En Bolivia, el fenómeno del pánico de compra se observó claramente cuando reportes en redes sociales sobre escasez generaron filas de horas y vaciaron surtidores en pocas horas.`
    : `La dinámica de vaciado es gradual y predecible. Este escenario corresponde a una operación con déficit moderado, similar al observado durante períodos de reducción de importaciones sin bloqueos activos.`;

  const metodosTexto = errEulerRK4 !== null
    ? `<strong>Comparativa de métodos:</strong> la diferencia entre Euler y RK4 al día ${dias} es <strong>${errEulerRK4} L</strong>. Esta diferencia se acumula paso a paso: Euler usa solo la pendiente al inicio del intervalo (error global O(h)), mientras que RK4 usa 4 evaluaciones intermedias (error global O(h⁴)). Para decisiones críticas — como cuándo activar la alerta de reserva baja o cuándo iniciar el racionamiento — RK4 provee la fecha más confiable.`
    : `Para comparar la acumulación de errores entre los 3 métodos, selecciona la opción <em>"Comparar los 3 métodos"</em> en el selector de método.`;

  document.getElementById('m5-interpretation').innerHTML = `
    <strong>Interpretación:</strong>
    Con entrada de <strong>${entrada} L/día</strong> y consumo base de <strong>${consumo} L/día</strong>
    (factor de pánico κ = ${panico}), la reserva ${idxVacio>=0 ? `<strong>se agota completamente el día ${diaVacio}</strong>` : `<strong>no se agota en los ${dias} días simulados</strong>`}.
    El nivel de alerta crítica (15% de la reserva inicial) se alcanza el <strong>día ${diaCrit}</strong>.
    <br/><br/>
    ${escenarioTexto}
    <br/><br/>
    ${metodosTexto}
  `;
  document.getElementById('m5-results').style.display = 'block';
}
