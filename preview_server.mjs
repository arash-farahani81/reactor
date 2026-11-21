import http from "node:http";
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const PORT = Number(process.env.PORT || 3000);
const HOST = "0.0.0.0";

const ALLOWED_FIGURES = [
  {
    id: "tcr_multiphysics_profiles.png",
    title: "TCR Multiphysics Axial Profiles",
    category: "Core Profiles",
    caption: "Axial temperature, conversion, and heat flux profiles across Co-Current, Counter-Current, and Conventional CMR."
  },
  {
    id: "nsga2_pareto_front.png",
    title: "NSGA-II Multi-Objective Pareto Frontier",
    category: "Optimization",
    caption: "Non-dominated trade-off between Methanol Yield, Benzene Production Rate, and Peak Tube Temperature."
  },
  {
    id: "comparison_coupled_vs_conventional.png",
    title: "Coupled TCR vs. Conventional CMR",
    category: "Comparison",
    caption: "Side-by-side thermal and yield comparison at locked equal syngas throughput (0.640 mol/s per tube)."
  },
  {
    id: "sensitivity_analysis.png",
    title: "Parametric Sensitivity Sweeps",
    category: "Sensitivity",
    caption: "Parametric sensitivity over tube inlet temperature, shell feed rate, and catalyst effectiveness factor."
  },
  {
    id: "reactor_results.png",
    title: "Co-Current Reactor State Profiles",
    category: "Core Profiles",
    caption: "Detailed species molar flows, temperatures, and pressure profiles along the 7.022 m reactor."
  },
  {
    id: "figures/fig1_thermal_profiles.png",
    title: "Figure 1 · Thermal Profiles",
    category: "Manuscript Figures",
    caption: "Publication figure of tube and shell temperature trajectories across operating modes."
  },
  {
    id: "figures/fig2_conversion_profiles.png",
    title: "Figure 2 · Conversion Profiles",
    category: "Manuscript Figures",
    caption: "Publication figure of CO2, CO, and Cyclohexane conversion trajectories."
  },
  {
    id: "figures/fig3_enthalpy_balance.png",
    title: "Figure 3 · Stream Enthalpy Balance",
    category: "Manuscript Figures",
    caption: "Independent stream enthalpy closure and heat transfer verification."
  },
  {
    id: "figures/fig4_pareto_frontier.png",
    title: "Figure 4 · Pareto Frontier",
    category: "Manuscript Figures",
    caption: "Publication figure of the multi-objective NSGA-II Pareto surface and knee point."
  }
];

const ALLOWED_DOCS = new Set([
  "simulation_summary.md",
  "config.yaml",
  "TCR_Methanol_Cyclohexane_Manuscript.md"
]);

const FIGURE_MAP = new Map(ALLOWED_FIGURES.map((f) => [f.id, f]));

function sendJson(res, status, data) {
  const body = JSON.stringify(data);
  res.writeHead(status, {
    "Content-Type": "application/json; charset=utf-8",
    "Cache-Control": "no-store"
  });
  res.end(body);
}

function sendText(res, status, contentType, text) {
  res.writeHead(status, {
    "Content-Type": `${contentType}; charset=utf-8`,
    "Cache-Control": "no-store"
  });
  res.end(text);
}

const INDEX_HTML = `<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>TCR Simulation Studio · Methanol Synthesis &amp; Cyclohexane Dehydrogenation</title>
  <meta name="description" content="Thermally Coupled Shell-and-Tube Reactor (TCR) rigorous multiphysics simulation, NSGA-II optimization, and Khademi benchmark viewer." />
  <link rel="preconnect" href="https://fonts.googleapis.com" />
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin />
  <link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600&family=IBM+Plex+Sans:wght@400;500;600;700&display=swap" rel="stylesheet" />
  <style>
    :root {
      --bg-canvas: #f8fafc;
      --bg-surface: #ffffff;
      --bg-subtle: #f1f5f9;
      --border-hairline: #e2e8f0;
      --border-strong: #cbd5e1;
      --text-primary: #0f172a;
      --text-secondary: #475569;
      --text-muted: #64748b;
      --accent: #2563eb;
      --accent-hover: #1d4ed8;
      --accent-soft: #eff6ff;
      --ok: #047857;
      --warn: #b45309;
      --danger: #be123c;
    }
    * { box-sizing: border-box; }
    html, body {
      margin: 0;
      padding: 0;
      background: var(--bg-canvas);
      color: var(--text-primary);
      font-family: "IBM Plex Sans", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      font-size: 15px;
      line-height: 1.5;
      -webkit-font-smoothing: antialiased;
    }
    .mono {
      font-family: "IBM Plex Mono", ui-monospace, SFMono-Regular, Menlo, monospace;
      font-variant-numeric: tabular-nums;
    }
    header {
      position: sticky;
      top: 0;
      z-index: 20;
      background: var(--bg-surface);
      border-bottom: 1px solid var(--border-hairline);
    }
    .topbar {
      max-width: 1360px;
      margin: 0 auto;
      padding: 0 24px;
      height: 56px;
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 32px;
    }
    .brand {
      font-size: 16px;
      font-weight: 700;
      letter-spacing: -0.01em;
      color: var(--text-primary);
      text-decoration: none;
      white-space: nowrap;
      flex-shrink: 0;
    }
    nav {
      display: flex;
      align-items: center;
      gap: 24px;
    }
    .nav-btn {
      background: none;
      border: none;
      padding: 6px 0;
      font-family: inherit;
      font-size: 14px;
      font-weight: 500;
      color: var(--text-secondary);
      cursor: pointer;
      white-space: nowrap;
      flex-shrink: 0;
      border-bottom: 2px solid transparent;
      transition: color 150ms ease, border-color 150ms ease;
    }
    .nav-btn:hover { color: var(--text-primary); }
    .nav-btn.active {
      color: var(--accent);
      border-bottom-color: var(--accent);
    }
    .nav-btn:focus-visible, .btn-primary:focus-visible, .seg-btn:focus-visible {
      outline: 2px solid var(--accent);
      outline-offset: 2px;
    }
    .btn-primary {
      background: var(--accent);
      color: #ffffff;
      border: 1px solid var(--accent);
      border-radius: 6px;
      padding: 8px 16px;
      font-family: inherit;
      font-size: 13px;
      font-weight: 600;
      cursor: pointer;
      white-space: nowrap;
      flex-shrink: 0;
      transition: background 150ms ease;
    }
    .btn-primary:hover { background: var(--accent-hover); }
    main {
      max-width: 1360px;
      margin: 0 auto;
      padding: 28px 24px 64px;
    }
    .page-header {
      margin-bottom: 24px;
      padding-bottom: 20px;
      border-bottom: 1px solid var(--border-hairline);
      display: flex;
      flex-wrap: wrap;
      align-items: flex-end;
      justify-content: space-between;
      gap: 16px;
    }
    h1 {
      margin: 0 0 6px;
      font-size: clamp(1.5rem, 2.4vw, 2rem);
      font-weight: 700;
      letter-spacing: -0.02em;
      text-wrap: balance;
    }
    .meta-line {
      font-size: 13px;
      color: var(--text-muted);
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      align-items: center;
    }
    .kpi-grid {
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      gap: 16px;
      margin-bottom: 24px;
    }
    .panel {
      background: var(--bg-surface);
      border: 1px solid var(--border-hairline);
      border-radius: 8px;
      padding: 18px 20px;
    }
    .kpi-label {
      font-size: 12px;
      font-weight: 500;
      color: var(--text-muted);
      margin-bottom: 8px;
    }
    .kpi-value-row {
      display: flex;
      align-items: baseline;
      gap: 6px;
      margin-bottom: 6px;
    }
    .kpi-value {
      font-size: 28px;
      font-weight: 600;
      color: var(--text-primary);
      line-height: 1.1;
    }
    .kpi-unit {
      font-size: 12px;
      color: var(--text-muted);
    }
    .kpi-sub {
      font-size: 12px;
      color: var(--text-secondary);
    }
    .section-title {
      margin: 0 0 4px;
      font-size: 16px;
      font-weight: 600;
      color: var(--text-primary);
    }
    .section-sub {
      margin: 0 0 16px;
      font-size: 13px;
      color: var(--text-muted);
    }
    .table-wrap {
      overflow-x: auto;
    }
    table {
      width: 100%;
      border-collapse: collapse;
      font-size: 13.5px;
    }
    th {
      text-align: left;
      padding: 10px 12px;
      font-size: 12px;
      font-weight: 600;
      color: var(--text-muted);
      border-bottom: 1px solid var(--border-strong);
      white-space: nowrap;
    }
    th.num, td.num { text-align: right; }
    td {
      padding: 11px 12px;
      border-bottom: 1px solid var(--border-hairline);
      color: var(--text-primary);
      white-space: nowrap;
    }
    tbody tr.selectable {
      cursor: pointer;
      transition: background 120ms ease;
    }
    tbody tr.selectable:hover { background: var(--bg-subtle); }
    tbody tr.selectable.selected { background: var(--accent-soft); }
    .split-grid {
      display: grid;
      grid-template-columns: 1.55fr 1fr;
      gap: 20px;
      margin-top: 20px;
      align-items: start;
    }
    .kv-list {
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 12px 20px;
      margin-top: 12px;
    }
    .kv-item {
      padding-bottom: 10px;
      border-bottom: 1px solid var(--border-hairline);
    }
    .kv-label {
      font-size: 12px;
      color: var(--text-muted);
      margin-bottom: 2px;
    }
    .kv-val {
      font-size: 14px;
      font-weight: 600;
      color: var(--text-primary);
    }
    .status-ok { color: var(--ok); font-weight: 600; }
    .status-warn { color: var(--warn); font-weight: 600; }
    .status-danger { color: var(--danger); font-weight: 600; }
    .figures-layout {
      display: grid;
      grid-template-columns: 310px 1fr;
      gap: 20px;
      align-items: start;
    }
    .fig-list {
      display: flex;
      flex-direction: column;
      gap: 6px;
    }
    .seg-btn {
      width: 100%;
      text-align: left;
      background: var(--bg-surface);
      border: 1px solid var(--border-hairline);
      border-radius: 6px;
      padding: 10px 12px;
      font-family: inherit;
      cursor: pointer;
      transition: background 120ms ease, border-color 120ms ease;
    }
    .seg-btn:hover { background: var(--bg-subtle); }
    .seg-btn.active {
      background: var(--accent-soft);
      border-color: var(--accent);
    }
    .seg-title {
      display: block;
      font-size: 13px;
      font-weight: 600;
      color: var(--text-primary);
      margin-bottom: 2px;
    }
    .seg-cat {
      display: block;
      font-size: 12px;
      color: var(--text-muted);
    }
    .figure-frame {
      background: #ffffff;
      border: 1px solid var(--border-hairline);
      border-radius: 8px;
      padding: 20px;
    }
    .figure-frame img {
      width: 100%;
      height: auto;
      display: block;
      border: 1px solid var(--border-hairline);
      border-radius: 4px;
      background: #ffffff;
    }
    .bench-grid {
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 20px;
      margin-bottom: 20px;
    }
    .chart-grid {
      display: grid;
      grid-template-columns: repeat(3, minmax(0, 1fr));
      gap: 16px;
    }
    .doc-tabs {
      display: flex;
      gap: 8px;
      margin-bottom: 16px;
      flex-wrap: wrap;
    }
    .doc-tab {
      background: var(--bg-subtle);
      border: 1px solid var(--border-hairline);
      border-radius: 6px;
      padding: 8px 14px;
      font-family: inherit;
      font-size: 13px;
      font-weight: 500;
      color: var(--text-secondary);
      cursor: pointer;
      white-space: nowrap;
    }
    .doc-tab.active {
      background: var(--accent);
      border-color: var(--accent);
      color: #ffffff;
    }
    pre.doc-viewer {
      margin: 0;
      padding: 18px 20px;
      background: var(--bg-surface);
      border: 1px solid var(--border-hairline);
      border-radius: 8px;
      font-family: "IBM Plex Mono", ui-monospace, monospace;
      font-size: 12.5px;
      line-height: 1.6;
      overflow-x: auto;
      white-space: pre-wrap;
      word-break: break-word;
      max-height: 70vh;
      overflow-y: auto;
      color: var(--text-primary);
    }
    .hidden { display: none !important; }
    .empty-note {
      padding: 24px;
      color: var(--text-muted);
      font-size: 14px;
    }
    @media (max-width: 1080px) {
      .kpi-grid { grid-template-columns: repeat(2, minmax(0, 1fr)); }
      .split-grid, .figures-layout, .bench-grid, .chart-grid { grid-template-columns: 1fr; }
    }
    @media (max-width: 640px) {
      .topbar { gap: 12px; padding: 0 14px; }
      nav { gap: 12px; overflow-x: auto; }
      .kpi-grid { grid-template-columns: 1fr; }
    }
  </style>
</head>
<body>
  <header>
    <div class="topbar">
      <a href="#overview" class="brand" id="brand-link">TCR Simulation Studio</a>
      <nav aria-label="Primary">
        <button type="button" class="nav-btn active" data-view="overview">Overview</button>
        <button type="button" class="nav-btn" data-view="figures">Figures</button>
        <button type="button" class="nav-btn" data-view="benchmark">Benchmark</button>
        <button type="button" class="nav-btn" data-view="files">Artifacts</button>
      </nav>
      <button type="button" class="btn-primary" id="refresh-btn">Refresh Data</button>
    </div>
  </header>

  <main>
    <div class="page-header">
      <div>
        <h1>Thermally Coupled Reactor Simulation Workbench</h1>
        <div class="meta-line" id="meta-line">
          <span>Methanol Synthesis + Cyclohexane Dehydrogenation</span>
          <span aria-hidden="true">&middot;</span>
          <span class="mono" id="meta-tubes">2962 tubes</span>
          <span aria-hidden="true">&middot;</span>
          <span class="mono" id="meta-updated">Loading results...</span>
        </div>
      </div>
    </div>

    <section id="view-overview">
      <div class="kpi-grid" id="kpi-grid"></div>
      <div class="panel">
        <h2 class="section-title">Reactor Configuration Comparison</h2>
        <p class="section-sub">Equal syngas throughput (0.640 mol/s per tube). Select any row to inspect conservation closures, thermal duties, and mechanical margins.</p>
        <div class="table-wrap">
          <table>
            <thead>
              <tr>
                <th>Configuration</th>
                <th class="num">MeOH Yield (%)</th>
                <th class="num">MeOH (t/d)</th>
                <th class="num">CH Conv. (%)</th>
                <th class="num">Benzene (t/d)</th>
                <th class="num">H2 (t/d)</th>
                <th class="num">T_max (K)</th>
                <th class="num">T_cl (K)</th>
                <th class="num">dP_shell (bar)</th>
                <th class="num">Steam Credit (MW)</th>
              </tr>
            </thead>
            <tbody id="regime-tbody"></tbody>
          </table>
        </div>
      </div>

      <div class="split-grid">
        <div class="panel">
          <h2 class="section-title" id="detail-title">Selected Regime Telemetry</h2>
          <p class="section-sub" id="detail-sub">Click a configuration row above to inspect detailed physical balances.</p>
          <div class="kv-list" id="detail-kv"></div>
        </div>
        <div class="panel">
          <h2 class="section-title">Plant Thermal Duties</h2>
          <p class="section-sub">Feed preheating, separation, and reaction-zone cooling utility summary.</p>
          <div class="kv-list" id="duties-kv"></div>
        </div>
      </div>
    </section>

    <section id="view-figures" class="hidden">
      <div class="figures-layout">
        <div class="fig-list" id="fig-list"></div>
        <div class="figure-frame">
          <h2 class="section-title" id="fig-title">Select a figure</h2>
          <p class="section-sub" id="fig-caption"></p>
          <img id="fig-img" alt="Simulation figure" referrerpolicy="no-referrer" />
        </div>
      </div>
    </section>

    <section id="view-benchmark" class="hidden">
      <div class="bench-grid">
        <div class="panel">
          <h2 class="section-title">Khademi et al. (2009) Benchmark Audit</h2>
          <p class="section-sub">Comparison against published OTCR/CMR benchmark and Table 4 equilibrium verification.</p>
          <div class="table-wrap">
            <table>
              <thead>
                <tr>
                  <th>Case</th>
                  <th class="num">MeOH Yield (%)</th>
                  <th class="num">CH Conv. (%)</th>
                  <th class="num">T_max (K)</th>
                </tr>
              </thead>
              <tbody id="bench-tbody"></tbody>
            </table>
          </div>
          <div class="kv-list" id="bench-notes"></div>
        </div>
        <div class="panel">
          <h2 class="section-title">OTP Ceiling &amp; EOS / Thiele Diagnostics</h2>
          <p class="section-sub">Theoretical optimal temperature profile ceiling, real-gas EOS delta, and Thiele modulus statistics.</p>
          <div class="kv-list" id="otp-kv"></div>
        </div>
      </div>

      <div class="panel">
        <h2 class="section-title">Parametric Sensitivity Curves</h2>
        <p class="section-sub">Rendered directly from <span class="mono">simulation_results.json</span> sensitivity sweeps.</p>
        <div class="chart-grid" id="sens-charts"></div>
      </div>
    </section>

    <section id="view-files" class="hidden">
      <div class="doc-tabs" id="doc-tabs">
        <button type="button" class="doc-tab active" data-doc="simulation_summary.md">simulation_summary.md</button>
        <button type="button" class="doc-tab" data-doc="config.yaml">config.yaml</button>
        <button type="button" class="doc-tab" data-doc="TCR_Methanol_Cyclohexane_Manuscript.md">TCR_Methanol_Cyclohexane_Manuscript.md</button>
      </div>
      <pre class="doc-viewer" id="doc-viewer">Loading document...</pre>
    </section>
  </main>

  <script>
    const state = {
      data: null,
      figures: [],
      selectedRegimeKey: "tcr_nsga2_optimum",
      selectedFigureId: null,
      selectedDoc: "simulation_summary.md",
      activeView: "overview"
    };

    const REGIMES = [
      { key: "conventional_baseline", label: "CMR Baseline" },
      { key: "conventional_optimized", label: "CMR Optimized" },
      { key: "cocurrent_baseline", label: "TCR Co-Current" },
      { key: "countercurrent_baseline", label: "TCR Counter-Current" },
      { key: "tcr_nsga2_optimum", label: "TCR NSGA-II Optimum" }
    ];

    function fmt(v, digits = 2, fallback = "—") {
      if (v === null || v === undefined || Number.isNaN(Number(v))) return fallback;
      return Number(v).toFixed(digits);
    }

    function fmtSci(v) {
      if (v === null || v === undefined || Number.isNaN(Number(v))) return "—";
      const n = Number(v);
      return n === 0 ? "0.00e+0" : n.toExponential(2);
    }

    function setView(viewName) {
      state.activeView = viewName;
      document.querySelectorAll(".nav-btn").forEach((b) => {
        b.classList.toggle("active", b.dataset.view === viewName);
      });
      ["overview", "figures", "benchmark", "files"].forEach((v) => {
        document.getElementById("view-" + v).classList.toggle("hidden", v !== viewName);
      });
      if (viewName === "files") loadDoc(state.selectedDoc);
    }

    document.querySelectorAll(".nav-btn").forEach((btn) => {
      btn.addEventListener("click", () => setView(btn.dataset.view));
    });
    document.getElementById("brand-link").addEventListener("click", (e) => {
      e.preventDefault();
      setView("overview");
    });
    document.getElementById("refresh-btn").addEventListener("click", () => loadAll());

    document.querySelectorAll(".doc-tab").forEach((btn) => {
      btn.addEventListener("click", () => {
        state.selectedDoc = btn.dataset.doc;
        document.querySelectorAll(".doc-tab").forEach((t) => {
          t.classList.toggle("active", t.dataset.doc === state.selectedDoc);
        });
        loadDoc(state.selectedDoc);
      });
    });

    async function loadDoc(name) {
      const el = document.getElementById("doc-viewer");
      el.textContent = "Loading " + name + "...";
      try {
        const resp = await fetch("/api/file?name=" + encodeURIComponent(name));
        if (!resp.ok) throw new Error("HTTP " + resp.status);
        el.textContent = await resp.text();
      } catch (err) {
        el.textContent = "Failed to load " + name + ": " + err.message;
      }
    }

    function renderOverview() {
      const d = state.data;
      if (!d) return;
      const meta = d.metadata || {};
      document.getElementById("meta-tubes").textContent = (meta.N_tubes || 2962) + " tubes · L = 7.022 m";
      document.getElementById("meta-updated").textContent =
        "SSOT v" + (meta.version || "2.0.0") + " · Generated " + (meta.generated_utc || "N/A");

      const opt = d.tcr_nsga2_optimum || {};
      const cmrOpt = d.conventional_optimized || {};
      const cmrBase = d.conventional_baseline || {};
      const steamCredit = cmrBase.cmr_steam_credit_MW ?? (cmrBase.Q_wall_total_MW ? Math.abs(cmrBase.Q_wall_total_MW) : null);

      const kpis = [
        {
          label: "NSGA-II Optimum MeOH Yield",
          value: fmt(opt.yield_net_COx ?? opt.yield_khademi_eq28, 2),
          unit: "%",
          sub: fmt(opt.daily_meoh_ton, 2) + " t/d vs. " + fmt(cmrOpt.daily_meoh_ton, 2) + " t/d CMR Opt"
        },
        {
          label: "Coproduced Benzene & H2",
          value: fmt(opt.daily_benzene_ton, 2),
          unit: "t/d Bz",
          sub: fmt(opt.daily_h2_ton, 2) + " t/d H2 (" + fmt(opt.H2_rate_plant, 2) + " mol/s)"
        },
        {
          label: "Peak Bed / Centerline Temp",
          value: fmt(opt.T_max_exo, 2),
          unit: "K",
          sub: "Centerline " + fmt(opt.T_centerline_max, 2) + " K · Limit 535.0 K"
        },
        {
          label: "CMR Steam Credit / TCR Cooling",
          value: fmt(steamCredit, 2),
          unit: "MW",
          sub: "TCR reaction-zone external cooling: 0.0 kW"
        }
      ];

      document.getElementById("kpi-grid").innerHTML = kpis.map((k) =>
        '<div class="panel">' +
          '<div class="kpi-label">' + k.label + '</div>' +
          '<div class="kpi-value-row">' +
            '<span class="kpi-value mono">' + k.value + '</span>' +
            '<span class="kpi-unit mono">' + k.unit + '</span>' +
          '</div>' +
          '<div class="kpi-sub">' + k.sub + '</div>' +
        '</div>'
      ).join("");

      const tbody = document.getElementById("regime-tbody");
      tbody.innerHTML = REGIMES.map((r) => {
        const m = d[r.key] || {};
        const isConv = r.key.startsWith("conventional");
        const sc = m.cmr_steam_credit_MW !== undefined
          ? m.cmr_steam_credit_MW
          : (isConv && m.Q_wall_total_MW !== undefined ? Math.abs(m.Q_wall_total_MW) : 0.0);
        const sel = state.selectedRegimeKey === r.key ? " selected" : "";
        return (
          '<tr class="selectable' + sel + '" data-key="' + r.key + '">' +
            '<td><strong>' + r.label + '</strong></td>' +
            '<td class="num mono">' + fmt(m.yield_net_COx ?? m.yield_khademi_eq28, 2) + '</td>' +
            '<td class="num mono">' + fmt(m.daily_meoh_ton, 2) + '</td>' +
            '<td class="num mono">' + (isConv ? "—" : fmt(m.X_CH, 2)) + '</td>' +
            '<td class="num mono">' + (isConv ? "0.00" : fmt(m.daily_benzene_ton, 2)) + '</td>' +
            '<td class="num mono">' + (isConv ? "0.00" : fmt(m.daily_h2_ton, 2)) + '</td>' +
            '<td class="num mono">' + fmt(m.T_max_exo, 2) + '</td>' +
            '<td class="num mono">' + fmt(m.T_centerline_max, 2) + '</td>' +
            '<td class="num mono">' + (isConv ? "—" : fmt(m.dP_shell, 3)) + '</td>' +
            '<td class="num mono">' + fmt(sc, 2) + '</td>' +
          '</tr>'
        );
      }).join("");

      tbody.querySelectorAll("tr[data-key]").forEach((row) => {
        row.addEventListener("click", () => {
          state.selectedRegimeKey = row.dataset.key;
          renderOverview();
        });
      });

      const activeMeta = REGIMES.find((r) => r.key === state.selectedRegimeKey) || REGIMES[4];
      const m = d[state.selectedRegimeKey] || {};
      document.getElementById("detail-title").textContent = activeMeta.label + " · Verification Telemetry";
      document.getElementById("detail-sub").textContent =
        "Mode: " + (m.mode || "N/A") + " · F_exo = " + fmt(m.F_exo_total_0 ?? 0.64, 3) + " mol/s per tube";

      const liftRatio = m.bed_lifting_ratio !== undefined
        ? m.bed_lifting_ratio
        : (m.dP_shell ? m.dP_shell / 0.74166 : 0.0);
      const needsGrid = m.requires_hold_down_grid !== undefined
        ? m.requires_hold_down_grid
        : liftRatio >= 0.85;
      const tClOk = (m.T_centerline_max || 0) <= 535.05;
      const elemC = m.elem_err_tube ? m.elem_err_tube.C : null;

      const detailItems = [
        ["CO2 / CO Conversion", fmt(m.X_CO2, 2) + "% / " + fmt(m.X_CO, 2) + "%"],
        ["Wall Heat Transfer Q_wall", fmt(m.Q_wall_total_MW, 3) + " MW"],
        ["Tube / Shell dP", fmt(m.dP_tube, 3) + " / " + fmt(m.dP_shell, 3) + " bar"],
        ["Shell Exit Pressure", m.P_shell_out == null ? "N/A (Boiling Water)" : fmt(m.P_shell_out, 3) + " bar"],
        ["Tube Enthalpy Error", fmt(m.enthalpy_err_tube_pct, 4) + "%"],
        ["Tube Carbon Balance Error", fmtSci(elemC)],
        ["ASME Hoop Stress / SF", fmt(m.hoop_stress_MPa, 2) + " MPa · SF " + fmt(m.asme_safety_factor, 2)],
        ["Centerline Sintering Limit", tClOk ? '<span class="status-ok">Nominal (&le; 535 K)</span>' : '<span class="status-warn">Exceeds 535 K (' + fmt(m.T_centerline_max, 1) + ' K)</span>'],
        ["Bed Lifting Ratio (dP / W)", m.dP_shell ? fmt(liftRatio, 2) + (needsGrid ? ' · <span class="status-warn">Hold-Down Grid Required</span>' : ' · <span class="status-ok">Self-Stable</span>') : "N/A"],
        ["Clamp Violations", String(m.clamp_violations ?? 0)]
      ];

      document.getElementById("detail-kv").innerHTML = detailItems.map(([k, v]) =>
        '<div class="kv-item"><div class="kv-label">' + k + '</div><div class="kv-val mono">' + v + '</div></div>'
      ).join("");

      const pt = d.process_thermal_duties || {};
      const dutyItems = [
        ["Syngas Preheat Duty", fmt(m.syngas_preheat_duty_MW ?? pt.syngas_preheat_duty_MW, 2) + " MW"],
        ["Shell Vaporization & Preheat", fmt(m.cyclohexane_vaporization_and_preheat_duty_MW ?? pt.cyclohexane_vaporization_and_preheat_duty_MW, 2) + " MW"],
        ["Total Feed Preheat", fmt(m.total_feed_preheat_duty_MW ?? pt.total_feed_preheat_duty_MW, 2) + " MW"],
        ["Extractive Distillation Reboiler", fmt(m.Q_reboiler_sep_MW, 2) + " MW"],
        ["H2 PSA Purification Work", fmt(m.W_psa_H2_MW, 2) + " MW"],
        ["Flash Condenser Total", fmt(m.Q_flash_condenser_MW, 2) + " MW"]
      ];
      document.getElementById("duties-kv").innerHTML = dutyItems.map(([k, v]) =>
        '<div class="kv-item"><div class="kv-label">' + k + '</div><div class="kv-val mono">' + v + '</div></div>'
      ).join("");
    }

    function renderFigures() {
      const listEl = document.getElementById("fig-list");
      if (!state.figures.length) {
        listEl.innerHTML = '<div class="empty-note">No PNG figures found in workspace.</div>';
        return;
      }
      if (!state.selectedFigureId) state.selectedFigureId = state.figures[0].id;

      listEl.innerHTML = state.figures.map((f) => {
        const active = f.id === state.selectedFigureId ? " active" : "";
        return (
          '<button type="button" class="seg-btn' + active + '" data-fig="' + f.id + '">' +
            '<span class="seg-title">' + f.title + '</span>' +
            '<span class="seg-cat">' + f.category + '</span>' +
          '</button>'
        );
      }).join("");

      listEl.querySelectorAll("button[data-fig]").forEach((btn) => {
        btn.addEventListener("click", () => {
          state.selectedFigureId = btn.dataset.fig;
          renderFigures();
        });
      });

      const cur = state.figures.find((f) => f.id === state.selectedFigureId) || state.figures[0];
      document.getElementById("fig-title").textContent = cur.title;
      document.getElementById("fig-caption").textContent = cur.caption;
      const img = document.getElementById("fig-img");
      img.src = "/img/" + cur.id + "?t=" + Date.now();
      img.alt = cur.title;
    }

    function makeMiniSvgChart(title, xLabel, xVals, seriesList) {
      if (!xVals || !xVals.length) {
        return '<div class="panel"><h3 class="section-title">' + title + '</h3><p class="section-sub">No sweep data available.</p></div>';
      }
      const W = 380, H = 210, padL = 48, padR = 16, padT = 16, padB = 34;
      const allY = seriesList.flatMap((s) => s.vals).filter((v) => typeof v === "number" && !Number.isNaN(v));
      const xMin = Math.min(...xVals), xMax = Math.max(...xVals);
      let yMin = Math.min(...allY), yMax = Math.max(...allY);
      if (Math.abs(yMax - yMin) < 1e-6) { yMin -= 1; yMax += 1; }
      const padY = (yMax - yMin) * 0.12;
      yMin -= padY; yMax += padY;

      const sx = (x) => padL + ((x - xMin) / Math.max(xMax - xMin, 1e-9)) * (W - padL - padR);
      const sy = (y) => H - padB - ((y - yMin) / Math.max(yMax - yMin, 1e-9)) * (H - padT - padB);

      const gridLines = [0, 0.25, 0.5, 0.75, 1].map((t) => {
        const yv = yMin + t * (yMax - yMin);
        const py = sy(yv).toFixed(1);
        return (
          '<line x1="' + padL + '" y1="' + py + '" x2="' + (W - padR) + '" y2="' + py + '" stroke="#e2e8f0" stroke-width="1" />' +
          '<text x="' + (padL - 6) + '" y="' + (Number(py) + 3.5) + '" text-anchor="end" font-size="10" fill="#64748b" font-family="IBM Plex Mono, monospace">' + yv.toFixed(1) + '</text>'
        );
      }).join("");

      const paths = seriesList.map((s) => {
        const pts = xVals.map((x, i) => sx(x).toFixed(1) + "," + sy(s.vals[i]).toFixed(1)).join(" ");
        const dots = xVals.map((x, i) =>
          '<circle cx="' + sx(x).toFixed(1) + '" cy="' + sy(s.vals[i]).toFixed(1) + '" r="3" fill="' + s.color + '" />'
        ).join("");
        return '<polyline fill="none" stroke="' + s.color + '" stroke-width="2" points="' + pts + '" />' + dots;
      }).join("");

      const legend = seriesList.map((s) =>
        '<span style="color:' + s.color + ';font-weight:600">' + s.label + '</span>'
      ).join(' <span style="color:#94a3b8">&middot;</span> ');

      return (
        '<div class="panel" style="padding:14px 16px">' +
          '<div style="display:flex;justify-content:space-between;align-items:baseline;gap:8px;margin-bottom:6px">' +
            '<span style="font-size:13px;font-weight:600">' + title + '</span>' +
            '<span style="font-size:11.5px">' + legend + '</span>' +
          '</div>' +
          '<svg viewBox="0 0 ' + W + ' ' + H + '" style="width:100%;height:auto;display:block" role="img" aria-label="' + title + '">' +
            gridLines +
            '<line x1="' + padL + '" y1="' + (H - padB) + '" x2="' + (W - padR) + '" y2="' + (H - padB) + '" stroke="#cbd5e1" stroke-width="1" />' +
            paths +
            '<text x="' + padL + '" y="' + (H - 12) + '" font-size="10" fill="#64748b" font-family="IBM Plex Mono, monospace">' + xMin + '</text>' +
            '<text x="' + ((padL + W - padR) / 2) + '" y="' + (H - 12) + '" text-anchor="middle" font-size="10.5" fill="#475569">' + xLabel + '</text>' +
            '<text x="' + (W - padR) + '" y="' + (H - 12) + '" text-anchor="end" font-size="10" fill="#64748b" font-family="IBM Plex Mono, monospace">' + xMax + '</text>' +
          '</svg>' +
        '</div>'
      );
    }

    function renderBenchmark() {
      const d = state.data;
      if (!d) return;
      const b = d.benchmark_khademi_2009 || {};
      const pub = b.khademi_published || {};
      const cmr = b.simulation_cmr || {};
      const thermo = b.simulation_tcr_thermo || {};
      const lit = b.simulation_tcr_literal || {};

      const rows = [
        ["Khademi (2009) Published OTCR", pub.tcr_yield_pct, pub.ch_conv_pct, pub.tcr_T_max_K],
        ["Simulated Benchmark CMR", cmr.yield_pct, null, cmr.T_max_K],
        ["Simulated OTCR (Thermo Kp)", thermo.yield_pct, thermo.X_CH_pct, thermo.T_max_K],
        ["Simulated OTCR (Literal Table 4)", lit.yield_pct, lit.X_CH_pct, lit.T_max_K]
      ];
      document.getElementById("bench-tbody").innerHTML = rows.map(([name, y, xch, tmax]) =>
        '<tr><td>' + name + '</td>' +
        '<td class="num mono">' + fmt(y, 2) + '</td>' +
        '<td class="num mono">' + fmt(xch, 2) + '</td>' +
        '<td class="num mono">' + fmt(tmax, 2) + '</td></tr>'
      ).join("");

      const typo = b.table4_typo_proof || {};
      document.getElementById("bench-notes").innerHTML = [
        ["Table 4 Printed -dH/R", fmt(typo.apparent_delta_H_over_R, 1) + " K"],
        ["Thermodynamic -dH/R", fmt(typo.true_delta_H_over_R, 1) + " K"]
      ].map(([k, v]) =>
        '<div class="kv-item"><div class="kv-label">' + k + '</div><div class="kv-val mono">' + v + '</div></div>'
      ).join("");

      const otp = b.otp_and_eos_audit || null;
      const th = b.thiele_analysis || {};
      const otpItems = [];
      if (otp) {
        const o035 = otp.otp_eta_035 || otp["eta_0.35"] || {};
        const o100 = otp.otp_eta_100 || otp["eta_1.0"] || {};
        const eos = otp.eos_comparison || otp.real_gas_vs_ideal || {};
        otpItems.push(["OTP Yield Ceiling (eta = 0.35)", fmt(o035.otp_yield_net, 2) + "% (Iso " + fmt(o035.best_isothermal_yield_net, 2) + "%)"]);
        otpItems.push(["OTP Yield Ceiling (eta = 1.00)", fmt(o100.otp_yield_net, 2) + "% (Iso " + fmt(o100.best_isothermal_yield_net, 2) + "%)"]);
        otpItems.push(["Real-Gas PR-EOS vs Ideal", fmt(eos.yield_net_real_eos, 2) + "% vs " + fmt(eos.yield_net_ideal, 2) + "%"]);
        otpItems.push(["EOS Yield Difference", fmt(eos.diff_pct_points, 3) + " pp"]);
      } else {
        otpItems.push(["OTP / EOS Audit", "Re-run run_all.py locally to populate"]);
      }
      otpItems.push(["Tube Thiele Mean / Min eta", fmt(th.tube_eta_mean, 3) + " / " + fmt(th.tube_eta_min, 3)]);
      otpItems.push(["Shell Thiele Mean / Min eta", fmt(th.shell_eta_mean, 3) + " / " + fmt(th.shell_eta_min, 3)]);
      document.getElementById("otp-kv").innerHTML = otpItems.map(([k, v]) =>
        '<div class="kv-item"><div class="kv-label">' + k + '</div><div class="kv-val mono">' + v + '</div></div>'
      ).join("");

      const sens = d.sensitivity_analysis || {};
      const tSw = sens.T_exo_sweep || {};
      const fSw = sens.F_endo_sweep || {};
      let etaX = [], etaCO2 = [], etaMeOH = [];
      if (Array.isArray(sens.eta_sweep)) {
        etaX = sens.eta_sweep.map((r) => r.eta);
        etaCO2 = sens.eta_sweep.map((r) => r.X_CO2);
        etaMeOH = sens.eta_sweep.map((r) => r.yield_meoh);
      } else if (sens.eta_sweep && sens.eta_sweep.eta_vals) {
        etaX = sens.eta_sweep.eta_vals;
        etaCO2 = sens.eta_sweep.X_CO2 || [];
        etaMeOH = sens.eta_sweep.yield_meoh || [];
      }

      document.getElementById("sens-charts").innerHTML = [
        makeMiniSvgChart("Tube Inlet Temp Sweep", "T_exo,0 (K)", tSw.T_vals, [
          { label: "X_CO2 (%)", vals: tSw.X_CO2 || [], color: "#2563eb" },
          { label: "X_CH (%)", vals: tSw.X_CH || [], color: "#0f766e" }
        ]),
        makeMiniSvgChart("Shell Feed Rate Sweep", "F_endo,0 (mol/s)", fSw.F_vals, [
          { label: "X_CO2 (%)", vals: fSw.X_CO2 || [], color: "#2563eb" },
          { label: "X_CH (%)", vals: fSw.X_CH || [], color: "#0f766e" }
        ]),
        makeMiniSvgChart("Catalyst Effectiveness Sweep", "eta_exo (-)", etaX, [
          { label: "MeOH Yield (%)", vals: etaMeOH, color: "#2563eb" },
          { label: "X_CO2 (%)", vals: etaCO2, color: "#b45309" }
        ])
      ].join("");
    }

    async function loadAll() {
      try {
        const [resResp, figResp] = await Promise.all([
          fetch("/api/results"),
          fetch("/api/figures")
        ]);
        state.data = await resResp.json();
        state.figures = await figResp.json();
        renderOverview();
        renderFigures();
        renderBenchmark();
      } catch (err) {
        document.getElementById("meta-updated").textContent = "Error loading simulation_results.json: " + err.message;
      }
    }

    loadAll();
  </script>
</body>
</html>`;

const server = http.createServer((req, res) => {
  const urlObj = new URL(req.url || "/", `http://${req.headers.host || "localhost"}`);
  const pathname = urlObj.pathname;

  if (pathname === "/" || pathname === "/index.html") {
    return sendText(res, 200, "text/html", INDEX_HTML);
  }

  if (pathname === "/api/results") {
    const target = path.join(__dirname, "simulation_results.json");
    try {
      const raw = fs.readFileSync(target, "utf8");
      return sendText(res, 200, "application/json", raw);
    } catch (err) {
      return sendJson(res, 500, { error: "Could not read simulation_results.json", detail: err.message });
    }
  }

  if (pathname === "/api/figures") {
    const available = ALLOWED_FIGURES.filter((f) => fs.existsSync(path.join(__dirname, f.id)));
    return sendJson(res, 200, available);
  }

  if (pathname === "/api/file") {
    const name = urlObj.searchParams.get("name") || "";
    if (!ALLOWED_DOCS.has(name)) {
      return sendJson(res, 404, { error: "Document not found" });
    }
    const target = path.join(__dirname, name);
    try {
      const raw = fs.readFileSync(target, "utf8");
      return sendText(res, 200, "text/plain", raw);
    } catch (err) {
      return sendJson(res, 500, { error: "Failed to read document", detail: err.message });
    }
  }

  if (pathname.startsWith("/img/")) {
    const figId = decodeURIComponent(pathname.slice("/img/".length));
    if (!FIGURE_MAP.has(figId)) {
      return sendJson(res, 404, { error: "Figure not found" });
    }
    const target = path.join(__dirname, figId);
    try {
      const buf = fs.readFileSync(target);
      res.writeHead(200, {
        "Content-Type": "image/png",
        "Cache-Control": "no-store"
      });
      res.end(buf);
      return;
    } catch (err) {
      return sendJson(res, 500, { error: "Failed to read figure", detail: err.message });
    }
  }

  sendJson(res, 404, { error: "Not found" });
});

server.listen(PORT, HOST, () => {
  console.log(`TCR Simulation Studio preview server running on http://${HOST}:${PORT}`);
});
