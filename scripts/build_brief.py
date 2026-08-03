# SPDX-License-Identifier: AGPL-3.0-or-later
"""Assemble the operational brief HTML (Artifact content + standalone for PDF)."""
import base64, json, os

OUT = os.path.join(os.path.dirname(__file__), "..", "docs", "fire_calvana_2026-07-03")
m = json.load(open(os.path.join(OUT, "metrics.json")))
mapb64 = base64.b64encode(open(os.path.join(OUT, "brief_map.png"), "rb").read()).decode()
atmob64 = base64.b64encode(open(os.path.join(OUT, "atmo_metrics_timeseries.png"), "rb").read()).decode()

iv = m["intervals"]
areas = [x["burned_ha"] for x in iv]
ros = [x["ros_head_p95_mmin"] for x in iv]
fli = [x["fli_max_kwm"] for x in iv]
peak_ros = max(ros); peak_fli = max(fli)
fli_typ = sorted(fli)[len(fli)//2]           # median cell-peak (typical)
conv = m["convection"]
ft = m["fire_type_cells"]
tot_cells = sum(int(v) for v in ft.values())
active_pct = round(100*int(ft.get("2",0))/tot_cells)
crown_pct = round(100*(int(ft.get("1",0))+int(ft.get("2",0)))/tot_cells)

# fire-heading text
head_line = next((l for l in m["operative"].splitlines() if "heading" in l), "")
op_rows = []
for name in ("head","right flank","tail","left flank"):
    for l in m["operative"].splitlines():
        s=l.strip()
        if s.startswith(name):
            import re
            mo=re.search(r"ROS\s+([\d.]+) ft/min.*driver:\s+(\w+)", s)
            if mo:
                op_rows.append((name, round(float(mo.group(1))*0.3048,1), mo.group(2)))
            break

def ember_bg(v, lo, hi):
    """map value -> rgba ember tint for intensity cells."""
    import colorsys
    t = max(0.0, min(1.0, (v-lo)/(hi-lo)))
    # amber(0)->deep red(1)
    stops = [(0.0,(242,166,59)),(0.5,(226,87,31)),(1.0,(161,20,20))]
    for (a,ca),(b,cb) in zip(stops, stops[1:]):
        if t<=b:
            f=(t-a)/(b-a) if b>a else 0
            r=int(ca[0]+(cb[0]-ca[0])*f); g=int(ca[1]+(cb[1]-ca[1])*f); bl=int(ca[2]+(cb[2]-ca[2])*f)
            return f"rgba({r},{g},{bl},0.22)"
    return "transparent"

rows_html = ""
for x in iv:
    rows_html += (
        f"<tr><td class='mono'>{x['clock_local']}</td>"
        f"<td class='mono num'>{x['burned_ha']:,.0f}</td>"
        f"<td class='mono num'>{x['ros_head_p95_mmin']:.0f}</td>"
        f"<td class='mono num' style='background:{ember_bg(x['fli_max_kwm'],40000,80000)}'>{x['fli_max_kwm']:,.0f}</td>"
        f"<td class='mono num'>{x['flame_max_m']:.0f}</td>"
        f"<td class='mono num'>{x['wind_ms']:.1f}</td>"
        f"<td class='mono num'>{x['rh_pct']:.0f}</td>"
        f"<td class='mono num'>{x['m1h_pct']:.0f}</td>"
        f"<td class='mono num'>{x['spots_cum'] if x['spots_cum'] is not None else '—'}</td></tr>\n")

op_html = "".join(
    f"<tr><td>{n}</td><td class='mono num'>{r}</td><td><span class='chip chip-wind'>{drv}</span></td></tr>"
    for n,r,drv in op_rows)

CSS = """
:root{
  --paper:#faf7f2; --raised:#fffdf9; --ink:#1c1712; --muted:#6f675d;
  --line:#e6ddd0; --ember:#e2571f; --amber:#f2a63b; --red:#a11414;
  --crit-bg:#f7e3dc; --crit-ink:#8a1a12; --warn-bg:#fbeecf; --warn-ink:#7a5410;
  --shadow:0 1px 2px rgba(60,40,20,.06),0 6px 20px rgba(60,40,20,.06);
}
@media (prefers-color-scheme:dark){:root{
  --paper:#17120d; --raised:#211a13; --ink:#ece3d6; --muted:#a89b8a;
  --line:#352a1f; --crit-bg:#3a1a14; --crit-ink:#f4b6a6; --warn-bg:#352915; --warn-ink:#f0cf90;
  --shadow:0 1px 2px rgba(0,0,0,.4),0 8px 24px rgba(0,0,0,.34);
}}
:root[data-theme="light"]{
  --paper:#faf7f2; --raised:#fffdf9; --ink:#1c1712; --muted:#6f675d; --line:#e6ddd0;
  --crit-bg:#f7e3dc; --crit-ink:#8a1a12; --warn-bg:#fbeecf; --warn-ink:#7a5410;
}
:root[data-theme="dark"]{
  --paper:#17120d; --raised:#211a13; --ink:#ece3d6; --muted:#a89b8a; --line:#352a1f;
  --crit-bg:#3a1a14; --crit-ink:#f4b6a6; --warn-bg:#352915; --warn-ink:#f0cf90;
}
*{box-sizing:border-box}
.wrap{max-width:1000px;margin:0 auto;padding:38px 26px 64px;
  background:var(--paper);color:var(--ink);
  font-family:system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;
  font-size:16px;line-height:1.55;-webkit-font-smoothing:antialiased;}
.serif{font-family:Georgia,"Iowan Old Style","Times New Roman",serif;}
.mono{font-family:ui-monospace,"SF Mono","JetBrains Mono",Menlo,monospace;
  font-variant-numeric:tabular-nums;}
.eyebrow{font-size:12px;letter-spacing:.14em;text-transform:uppercase;color:var(--muted);
  font-weight:600;margin:0 0 6px;}
h1{font-size:34px;line-height:1.1;margin:0;font-weight:700;text-wrap:balance;}
.sub{color:var(--muted);margin:8px 0 0;font-size:15px;}
.stripe{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin:18px 0 8px;}
.badge{display:inline-flex;align-items:center;gap:7px;font-weight:700;font-size:13px;
  letter-spacing:.05em;padding:6px 12px;border-radius:999px;background:var(--crit-bg);color:var(--crit-ink);}
.badge::before{content:"";width:8px;height:8px;border-radius:50%;background:var(--red);}
.meta{color:var(--muted);font-size:13px;}
.rule{height:3px;background:linear-gradient(90deg,var(--amber),var(--ember),var(--red));
  border:0;border-radius:3px;margin:20px 0 26px;}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:14px;margin:0 0 30px;}
.tile{background:var(--raised);border:1px solid var(--line);border-radius:12px;padding:16px 16px 14px;box-shadow:var(--shadow);}
.tile .k{font-size:12px;letter-spacing:.06em;text-transform:uppercase;color:var(--muted);font-weight:600;}
.tile .v{font-size:29px;font-weight:700;margin-top:6px;letter-spacing:-.01em;}
.tile .u{font-size:13px;color:var(--muted);font-weight:500;margin-left:3px;}
.tile .n{font-size:12px;color:var(--muted);margin-top:2px;}
h2{font-size:13px;letter-spacing:.12em;text-transform:uppercase;color:var(--ember);
  font-weight:700;margin:38px 0 14px;padding-bottom:8px;border-bottom:1px solid var(--line);}
figure{margin:0 0 8px;}
figure img{width:100%;height:auto;border-radius:12px;border:1px solid var(--line);display:block;box-shadow:var(--shadow);}
figcaption{font-size:13px;color:var(--muted);margin-top:9px;}
.scroll{overflow-x:auto;border:1px solid var(--line);border-radius:12px;box-shadow:var(--shadow);}
table{border-collapse:collapse;width:100%;font-size:14px;background:var(--raised);min-width:640px;}
th,td{padding:9px 12px;text-align:left;border-bottom:1px solid var(--line);}
thead th{font-size:11px;letter-spacing:.05em;text-transform:uppercase;color:var(--muted);
  font-weight:600;background:var(--paper);position:sticky;top:0;}
td.num,th.num{text-align:right;}
tbody tr:last-child td{border-bottom:0;}
tbody tr:last-child{font-weight:700;}
.cols{display:grid;grid-template-columns:1fr 1fr;gap:22px;}
@media (max-width:680px){.cols{grid-template-columns:1fr;}h1{font-size:27px;}}
.card{background:var(--raised);border:1px solid var(--line);border-radius:12px;padding:18px 20px;box-shadow:var(--shadow);}
.card h3{margin:0 0 10px;font-size:15px;}
.kv{display:flex;justify-content:space-between;gap:12px;padding:5px 0;border-bottom:1px dashed var(--line);font-size:14px;}
.kv:last-child{border-bottom:0;}.kv .lbl{color:var(--muted);}
.chip{display:inline-block;font-size:11px;font-weight:700;letter-spacing:.03em;padding:2px 9px;border-radius:999px;}
.chip-wind{background:var(--warn-bg);color:var(--warn-ink);}
.callout{background:var(--warn-bg);color:var(--warn-ink);border-left:4px solid var(--amber);
  border-radius:0 10px 10px 0;padding:16px 20px;margin:8px 0 0;}
.callout h3{margin:0 0 8px;font-size:15px;color:var(--warn-ink);}
.callout ol{margin:0;padding-left:20px;}.callout li{margin:7px 0;font-size:14.5px;line-height:1.5;}
.prov{font-size:12.5px;color:var(--muted);line-height:1.7;margin-top:10px;}
.prov b{color:var(--ink);font-weight:600;}
footer{margin-top:34px;padding-top:16px;border-top:1px solid var(--line);font-size:12px;color:var(--muted);}
"""

BODY = f"""
<div class="wrap">
  <p class="eyebrow">pyflam · operational propagation brief</p>
  <h1 class="serif">Calvana ridge wildfire</h1>
  <p class="sub">Ignition 43.93561 N, 11.09635 E · above Vaiano / Prato, Tuscany · analysis window 14:00–20:00 CEST, 3 July 2026</p>
  <div class="stripe">
    <span class="badge">EXTREME · active crown fire</span>
    <span class="meta mono">heading 222° (SW) · wind-driven · 6-hour projection</span>
  </div>
  <hr class="rule"/>

  <div class="tiles">
    <div class="tile"><div class="k">Burned area · 6 h</div><div class="v mono">6,631<span class="u">ha</span></div><div class="n">from a point ignition</div></div>
    <div class="tile"><div class="k">Head run</div><div class="v mono">12.8<span class="u">km</span></div><div class="n">to SW, along the ridge</div></div>
    <div class="tile"><div class="k">Peak head ROS</div><div class="v mono">{peak_ros:.0f}<span class="u">m/min</span></div><div class="n">sustained 40–52 m/min</div></div>
    <div class="tile"><div class="k">Active crown</div><div class="v mono">{active_pct}<span class="u">%</span></div><div class="n">{crown_pct}% crown-involved</div></div>
    <div class="tile"><div class="k">Flame length</div><div class="v mono">12–16<span class="u">m</span></div><div class="n">head, actively crowning</div></div>
    <div class="tile"><div class="k">Fireline intensity</div><div class="v mono">~75k<span class="u">kW/m</span></div><div class="n">upper-bound · see caveats</div></div>
  </div>

  <h2>Fire-front progression</h2>
  <figure>
    <img alt="30-minute fire-front isochrones over shaded terrain, spreading southwest from the ignition" src="data:image/png;base64,{mapb64}"/>
    <figcaption>Time-of-arrival isochrones at 30-minute intervals (14:30 → 20:00 CEST). Colour ramps from amber (early) to deep red (late). The front runs west–southwest off the Calvana ridge under the north-easterly wind, with a southern spur down-terrain. Grid 100 m; ✶ = ignition.</figcaption>
  </figure>

  <h2>Behaviour every 30 minutes</h2>
  <div class="scroll">
  <table>
    <thead><tr>
      <th>local</th><th class="num">area (ha)</th><th class="num">head ROS<br>(m/min)</th>
      <th class="num">FLI max<br>(kW/m)</th><th class="num">flame (m)</th><th class="num">wind (m/s)</th>
      <th class="num">RH %</th><th class="num">1-h FM %</th><th class="num">spot fires</th>
    </tr></thead>
    <tbody>{rows_html}</tbody>
  </table>
  </div>
  <p class="meta" style="margin-top:9px;font-size:12.5px;">Head ROS is the 95th-percentile heading spread rate over cells burned in each step. FLI max is the hottest cell that step (Byram intensity). Weather is the ICON-2i column at the ignition; fuel moisture is the equilibrium dead 1-h moisture derived from it.</p>

  <div class="cols" style="margin-top:30px;">
    <div>
      <h2 style="margin-top:0;">Convective / pyro behaviour</h2>
      <div class="card">
        <div class="kv"><span class="lbl">Continuous-Haines</span><span class="mono">{conv['continuous_haines']} / 13</span></div>
        <div class="kv"><span class="lbl">Inverted-V sounding</span><span class="mono">{'yes' if conv['inverted_v'] else 'no'}</span></div>
        <div class="kv"><span class="lbl">Sfc dewpoint depression</span><span class="mono">{conv['sfc_dewpoint_depression_C']} °C</span></div>
        <div class="kv"><span class="lbl">Mid-level RH</span><span class="mono">{conv['mid_rh_pct']:.0f} %</span></div>
        <div class="kv"><span class="lbl">Lifting condensation level</span><span class="mono">{conv['lcl_m']:.0f} m</span></div>
        <div class="kv"><span class="lbl">Plume loft factor</span><span class="mono">×{conv['plume_factor_max']}</span></div>
        <div class="kv"><span class="lbl">pyroCb potential</span><span class="mono">unfavourable</span></div>
      </div>
      <p class="meta" style="font-size:12.5px;margin-top:10px;">A dry, moderately unstable afternoon: deep dry sub-cloud layer (high LCL, 19 °C dewpoint depression) but moist enough aloft that a plume-dominated pyroconvective column is <b>not</b> favoured. Expect vigorous convection columns over the head, short-range spotting — not a pyroCb day.</p>
    </div>
    <div>
      <h2 style="margin-top:0;">What is driving the front</h2>
      <div class="scroll">
      <table>
        <thead><tr><th>sector</th><th class="num">ROS (m/min)</th><th>dominant driver</th></tr></thead>
        <tbody>{op_html}</tbody>
      </table>
      </div>
      <p class="meta" style="font-size:12.5px;margin-top:10px;">Operative decomposition at the final perimeter (18:00 UTC). <b>Wind is the dominant driver in every sector</b> — the fire is wind-aligned toward the SW; slope and fuel-load gradients are secondary. As the afternoon wind eases (5.0 → 3.7 m/s) the head rate relaxes but crowning persists in the dense canopy.</p>
    </div>
  </div>

  <h2>Atmospheric drivers — 6 h before to 6 h after</h2>
  <figure>
    <img alt="Five stacked time-series panels of ICON-2i atmospheric drivers at the ignition from 08:00 to 02:00 CEST: pyroconvection indices, convective moisture, dead fuel moisture, 10 m wind, and 2 m relative humidity and temperature, with the simulation window shaded" src="data:image/png;base64,{atmob64}"/>
    <figcaption>ICON-2i column at the ignition, 08:00 CEST (3 Jul) → 02:00 CEST (4 Jul); the 14:00–20:00 CEST simulation window is shaded. Note the diurnal signal: pyroconvective potential (Continuous-Haines, plume loft) actually <b>peaks mid-morning and eases into the burn window</b>; the afternoon is fuel- and wind-driven — 1-h dead fuel moisture bottoms near 5.7% around 17:00 while the north-easterly wind backs and eases from ~8 to ~4 m/s.</figcaption>
  </figure>

  <h2>Reading this brief — two caveats</h2>
  <div class="callout">
    <h3>⚠ Interpret intensity as an upper bound; this is a decision-support projection</h3>
    <ol>
      <li><b>Fireline intensities (~55,000–75,000 kW/m, single-cell peaks higher) are an upper bound.</b> They are driven hard by the canopy bulk density and load, which here are <b>derived from GEDI canopy height + cover with fire-science heuristics</b> — modelled estimates, not field measurements. Head rate of spread was independently checked against Cruz-2005 and is physically consistent.</li>
      <li><b>This run uses the ICON-2i single-column wind, not the coupled CFD fire-plume.</b> Cruz-2005 crown fire, ICON-2i weather, ember spotting and the pyroconvection sounding are all included; the buoyant-plume feedback on local wind and lofting is not. The full OpenFOAM plume coupling is available and can be run on request.</li>
    </ol>
  </div>

  <h2>Methods &amp; data lineage</h2>
  <p class="prov">
    <b>Engine</b> pyflam 0.2.0.dev — minimum-travel-time front tracking on a Rothermel surface kernel with Cruz-2005 active-crown-fire coupling; Byram fireline intensity &amp; flame length; firebrand spotting from crown intensity.<br/>
    <b>Weather</b> ICON-2i 2.2 km convection-permitting model (MISTRAL / AgenziaItaliaMeteo open archive), run 2026-07-03 00Z; 10 m wind, 2 m T/dewpoint → equilibrium dead fuel moisture; 850/700/500 hPa T + RH → Continuous-Haines, inverted-V, LCL.<br/>
    <b>Fuels</b> Tuscany Scott &amp; Burgan 40 fuel model grid (ignition cell TU5 #165); slope/aspect/DEM at 100 m.<br/>
    <b>Canopy</b> GEDI-calibrated canopy height (Meta/WRI) + cover; canopy base height &amp; bulk density derived by heuristic (Cruz, Alexander &amp; Wakimoto 2003).<br/>
    <b>Domain</b> 100 m grid, 22 × 20 km, tightened to the modelled footprint · start 14:00 CEST (12:00 UTC) · 360 min · 30-min steps.
  </p>

  <footer>
    Generated by pyflam · propagation module · not an official incident forecast. Coordinates WGS-84; grid EPSG:3035.
  </footer>
</div>
"""

content = f"<style>{CSS}</style>\n{BODY}"
with open(os.path.join(OUT, "operational_brief.html"), "w") as fh:
    fh.write(content)

standalone = ("<!doctype html><html><head><meta charset='utf-8'>"
              "<meta name='viewport' content='width=device-width,initial-scale=1'>"
              "<title>Calvana wildfire — operational brief</title>"
              "<style>html,body{margin:0;padding:0;background:#faf7f2}</style>"
              f"</head><body>{content}</body></html>")
with open(os.path.join(OUT, "operational_brief_standalone.html"), "w") as fh:
    fh.write(standalone)
print("wrote operational_brief.html and operational_brief_standalone.html")
print("map bytes b64:", len(mapb64))
