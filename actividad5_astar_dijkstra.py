"""
ACTIVIDAD 5 (Sección 2.5 Búsqueda con heurísticas: A* / Dijkstra) - PROYECTO EN EQUIPO
Fundamentos de Inteligencia Artificial - Universidad Politécnica de Atlautla

EQUIPO
------
  * Suarez Alvarado Ehud Jockzan
  * Rosales Huesca David Fernando
  * Palacios Martinez Josafad
  * Sanchez Valero Crhystian Samuel

CASO REAL: ruteo GPS entre ciudades de México
---------------------------------------------
Dada una red de 54 ciudades (coordenadas GPS aproximadas), encontrar la ruta de menor
distancia entre un origen y un destino, como lo haría un navegador GPS.

MODELADO COMO PROBLEMA DE BÚSQUEDA
  Estado    : ciudad actual.
  Acciones  : viajar a una ciudad conectada por carretera.
  Coste     : distancia del tramo en km.
  Meta      : ciudad destino.
  Grafo     : no dirigido y ponderado. Cada ciudad se conecta con sus 4 ciudades más
              cercanas (k-vecinos) y el peso es  Haversine(a, b) x FACTOR_SINUOSIDAD.

SUPUESTOS (declarar en el informe): la red de carreteras es una APROXIMACIÓN
construida por proximidad, no el mapa oficial; FACTOR_SINUOSIDAD = 1.2 es un
parámetro supuesto para que la distancia por carretera sea mayor que la línea recta.
Las coordenadas son aproximadas y deben verificarse. Con datos reales de OpenStreetMap
(p. ej. osmnx) el mismo código de búsqueda funciona sin cambios.

ALGORITMOS
  Dijkstra : expande por menor coste acumulado g(n). Óptimo, sin información de la meta.
  A*       : expande por f(n) = g(n) + h(n). h(n) = distancia en línea recta (Haversine)
             hasta el destino. Si h es ADMISIBLE (nunca sobreestima) A* es óptimo.
  Dijkstra equivale a A* con h(n) = 0.

ADMISIBILIDAD: h(n) <= h*(n) (coste real mínimo). Aquí se cumple porque cada arista
pesa al menos su distancia en línea recta y, por la desigualdad triangular, ningún
camino puede ser más corto que la línea recta. El programa lo VERIFICA numéricamente
(admisibilidad y consistencia) y muestra qué pasa con heurísticas NO admisibles.

Librerías: numpy, matplotlib, networkx (pip) + heapq, math, time, csv, itertools,
           dataclasses, pathlib, typing (biblioteca estándar).
"""

from __future__ import annotations

import csv
import heapq
import itertools
import math
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec

EQUIPO = ["Suarez Alvarado Ehud Jockzan", "Rosales Huesca David Fernando",
          "Palacios Martinez Josafad", "Sanchez Valero Crhystian Samuel"]

# ----------------------------------------------------------------------------
# 1. DATOS: ciudades (lat, lon) aproximadas, en grados
# ----------------------------------------------------------------------------
CIUDADES = {
    "CDMX": (19.4326, -99.1332), "Guadalajara": (20.6597, -103.3496),
    "Monterrey": (25.6866, -100.3161), "Puebla": (19.0414, -98.2063),
    "Querétaro": (20.5888, -100.3899), "León": (21.1250, -101.6860),
    "Mérida": (20.9674, -89.5926), "Veracruz": (19.1738, -96.1342),
    "Tijuana": (32.5149, -117.0382), "Chihuahua": (28.6330, -106.0691),
    "Oaxaca": (17.0732, -96.7266), "Cancún": (21.1619, -86.8515),
    "Hermosillo": (29.0729, -110.9559), "Culiacán": (24.8091, -107.3940),
    "Torreón": (25.5428, -103.4068), "San Luis Potosí": (22.1565, -100.9855),
    "Aguascalientes": (21.8853, -102.2916), "Morelia": (19.7060, -101.1950),
    "Toluca": (19.2826, -99.6557), "Cuernavaca": (18.9242, -99.2216),
    "Villahermosa": (17.9892, -92.9475), "Tuxtla Gutiérrez": (16.7528, -93.1152),
    "Acapulco": (16.8531, -99.8237), "Mazatlán": (23.2494, -106.4111),
    "Saltillo": (25.4232, -101.0053), "Tampico": (22.2331, -97.8611),
    "Pachuca": (20.1011, -98.7591), "Campeche": (19.8301, -90.5349),
    "Durango": (24.0277, -104.6532), "Mexicali": (32.6245, -115.4523),
    "Ciudad Juárez": (31.6904, -106.4245), "Nuevo Laredo": (27.4763, -99.5164),
    "Reynosa": (26.0508, -98.2979), "Matamoros": (25.8697, -97.5027),
    "Ciudad Victoria": (23.7369, -99.1411), "Zacatecas": (22.7709, -102.5832),
    "Tepic": (21.5042, -104.8946), "Colima": (19.2433, -103.7250),
    "Chilpancingo": (17.5506, -99.5058), "Xalapa": (19.5438, -96.9102),
    "Tlaxcala": (19.3139, -98.2404), "Guanajuato": (21.0190, -101.2574),
    "Chetumal": (18.5001, -88.2961), "Ensenada": (31.8667, -116.5964),
    "Los Mochis": (25.7935, -108.9859), "Ciudad Obregón": (27.4828, -109.9304),
    "Nogales": (31.3086, -110.9422), "Piedras Negras": (28.7000, -100.5235),
    "Tapachula": (14.9050, -92.2646), "Ciudad Valles": (21.9953, -99.0140),
    "Poza Rica": (20.5331, -97.4592), "Coatzacoalcos": (18.1500, -94.4333),
    "Irapuato": (20.6767, -101.3563), "Celaya": (20.5230, -100.8157),
}
K_VECINOS = 4
FACTOR_SINUOSIDAD = 1.2
CASOS = [("Tijuana", "Cancún"), ("CDMX", "Monterrey"), ("Guadalajara", "Mérida"),
         ("Tapachula", "Nogales"), ("Puebla", "Querétaro")]
CASO_PRINCIPAL = ("Tijuana", "Cancún")


def haversine(p, q, radio_km: float = 6371.0) -> float:
    """Distancia de círculo máximo (km) entre dos puntos (lat, lon)."""
    lat1, lon1, lat2, lon2 = map(math.radians, (*p, *q))
    a = (math.sin((lat2 - lat1) / 2) ** 2
         + math.cos(lat1) * math.cos(lat2) * math.sin((lon2 - lon1) / 2) ** 2)
    return 2 * radio_km * math.asin(math.sqrt(a))


# ----------------------------------------------------------------------------
# 2. GRAFO DEL PROBLEMA
# ----------------------------------------------------------------------------
def construir_grafo() -> nx.Graph:
    """Conecta cada ciudad con sus K más cercanas; peso = Haversine x sinuosidad."""
    G = nx.Graph()
    G.add_nodes_from(CIUDADES)
    for a, pa in CIUDADES.items():
        cercanas = sorted((b for b in CIUDADES if b != a),
                          key=lambda b: haversine(pa, CIUDADES[b]))[:K_VECINOS]
        for b in cercanas:
            G.add_edge(a, b, weight=haversine(pa, CIUDADES[b]) * FACTOR_SINUOSIDAD)
    while not nx.is_connected(G):                  # garantizar conectividad
        comps = list(nx.connected_components(G))
        a, b = min(((x, y) for x in comps[0] for c in comps[1:] for y in c),
                   key=lambda par: haversine(CIUDADES[par[0]], CIUDADES[par[1]]))
        G.add_edge(a, b, weight=haversine(CIUDADES[a], CIUDADES[b]) * FACTOR_SINUOSIDAD)
    return G


def lista_adyacencia(G: nx.Graph) -> Dict[str, List[Tuple[str, float]]]:
    return {u: [(v, G[u][v]["weight"]) for v in G[u]] for u in G.nodes}


# ----------------------------------------------------------------------------
# 3. HEURÍSTICAS
# ----------------------------------------------------------------------------
Heuristica = Callable[[str, str], float]


def h_cero(n: str, meta: str) -> float:
    """h = 0: convierte A* en Dijkstra."""
    return 0.0


def h_recta(n: str, meta: str) -> float:
    """Distancia en línea recta (Haversine). ADMISIBLE y CONSISTENTE."""
    return haversine(CIUDADES[n], CIUDADES[meta])


def h_recta_sinuosa(n: str, meta: str) -> float:
    """Línea recta x factor de sinuosidad: más informada, también admisible aquí."""
    return FACTOR_SINUOSIDAD * haversine(CIUDADES[n], CIUDADES[meta])


def h_ponderada(w: float, base: Heuristica = h_recta_sinuosa) -> Heuristica:
    """Weighted A*: w·h. Para w > 1 deja de ser admisible (pierde optimalidad)."""
    return lambda n, meta: w * base(n, meta)


# ----------------------------------------------------------------------------
# 4. ALGORITMOS (cola de prioridad con heapq)
# ----------------------------------------------------------------------------
@dataclass
class Resultado:
    nombre: str
    camino: Optional[List[str]] = None
    costo: float = math.inf
    expandidos: int = 0
    generados: int = 0
    max_abierta: int = 0
    tiempo_s: float = 0.0
    orden: List[str] = field(default_factory=list)       # nodos expandidos (para graficar)


def _mejor_primero(adj, inicio: str, meta: str, h: Heuristica, nombre: str) -> Resultado:
    """
    Búsqueda de mejor-primero con f(n) = g(n) + h(n).
      h = 0  -> Dijkstra      |     h admisible -> A* óptimo.
    Admite reabrir nodos si se encuentra un g menor (necesario si h no es consistente).
    """
    t0 = time.perf_counter()
    res = Resultado(nombre)
    g = {inicio: 0.0}
    padre: Dict[str, Optional[str]] = {inicio: None}
    contador = itertools.count()                          # desempate estable en el heap
    cache_h: Dict[str, float] = {}                        # h(n) se calcula una sola vez por nodo

    def h_cacheada(n: str) -> float:
        if n not in cache_h:
            cache_h[n] = h(n, meta)
        return cache_h[n]

    abierta = [(h_cacheada(inicio), next(contador), 0.0, inicio)]
    res.generados = 1
    while abierta:
        res.max_abierta = max(res.max_abierta, len(abierta))
        f, _, g_u, u = heapq.heappop(abierta)
        if g_u > g[u] + 1e-12:                            # entrada obsoleta
            continue
        res.expandidos += 1
        res.orden.append(u)
        if u == meta:
            camino, n = [], meta
            while n is not None:
                camino.append(n); n = padre[n]
            res.camino, res.costo = camino[::-1], g[meta]
            break
        for v, w in adj[u]:
            ng = g_u + w
            if ng < g.get(v, math.inf) - 1e-12:
                g[v], padre[v] = ng, u
                heapq.heappush(abierta, (ng + h_cacheada(v), next(contador), ng, v))
                res.generados += 1
    res.tiempo_s = time.perf_counter() - t0
    return res


def dijkstra(adj, inicio, meta) -> Resultado:
    return _mejor_primero(adj, inicio, meta, h_cero, "Dijkstra")


def a_estrella(adj, inicio, meta, h: Heuristica = h_recta, nombre: str = "A*") -> Resultado:
    return _mejor_primero(adj, inicio, meta, h, nombre)


def medir(f, repeticiones: int = 20) -> Resultado:
    """Ejecuta varias veces y conserva el menor tiempo (reduce el ruido)."""
    mejor = None
    for _ in range(repeticiones):
        r = f()
        if mejor is None or r.tiempo_s < mejor.tiempo_s:
            mejor = r
    return mejor


# ----------------------------------------------------------------------------
# 5. VERIFICACIÓN DE ADMISIBILIDAD Y CONSISTENCIA
# ----------------------------------------------------------------------------
def verificar_heuristica(G: nx.Graph, h: Heuristica, tol: float = 1e-9):
    """
    Para TODOS los destinos posibles:
      admisible : h(n, meta) <= h*(n, meta)    (h* = coste real, por Dijkstra)
      consistente: h(u, meta) <= w(u,v) + h(v, meta)  para toda arista (u, v)
    Devuelve (viol_admisibilidad, viol_consistencia, razones h/h*, total pares).
    """
    viol_a = viol_c = total = 0
    razones = []
    for meta in G.nodes:
        h_real = nx.single_source_dijkstra_path_length(G, meta, weight="weight")
        for n, real in h_real.items():
            if n == meta:
                continue
            total += 1
            hv = h(n, meta)
            razones.append(hv / real)
            if hv > real + tol:
                viol_a += 1
        for u, v, d in G.edges(data=True):
            w = d["weight"]
            if h(u, meta) > w + h(v, meta) + tol or h(v, meta) > w + h(u, meta) + tol:
                viol_c += 1
    return viol_a, viol_c, np.array(razones), total


# ----------------------------------------------------------------------------
# 6. EXPERIMENTOS
# ----------------------------------------------------------------------------
def estudio_casos(adj):
    filas = []
    for o, d in CASOS:
        rd = medir(lambda: dijkstra(adj, o, d))
        ra = medir(lambda: a_estrella(adj, o, d))
        assert abs(rd.costo - ra.costo) < 1e-6, "A* debe igualar el óptimo de Dijkstra"
        filas.append((o, d, rd, ra))
    return filas


def todos_los_pares(G, adj, pesos_w=(1.0, 1.25, 1.5, 2.0, 3.0)):
    """Ejecuta los algoritmos en los N(N-1)/2 pares y verifica contra NetworkX."""
    datos = []
    barrido = {w: {"exp": [], "sub": [], "opt": 0} for w in pesos_w}
    for o, d in itertools.combinations(G.nodes, 2):
        ref = nx.dijkstra_path_length(G, o, d, weight="weight")          # verificación externa
        rd = dijkstra(adj, o, d)
        ra = a_estrella(adj, o, d, h_recta, "A* (recta)")
        rs = a_estrella(adj, o, d, h_recta_sinuosa, "A* (recta·1.2)")
        assert abs(rd.costo - ref) < 1e-6 and abs(ra.costo - ref) < 1e-6 and abs(rs.costo - ref) < 1e-6
        datos.append((o, d, haversine(CIUDADES[o], CIUDADES[d]), ref, rd, ra, rs))
        for w in pesos_w:
            rw = a_estrella(adj, o, d, h_ponderada(w))
            barrido[w]["exp"].append(rw.expandidos)
            barrido[w]["sub"].append(100 * (rw.costo / ref - 1))
            barrido[w]["opt"] += abs(rw.costo - ref) < 1e-6
    return datos, barrido


# ----------------------------------------------------------------------------
# 7. DASHBOARD
# ----------------------------------------------------------------------------
def _segmentos(G):
    return [[(CIUDADES[u][1], CIUDADES[u][0]), (CIUDADES[v][1], CIUDADES[v][0])] for u, v in G.edges]


def dibujar_mapa(ax, G, res: Resultado, color_expandidos, titulo, o, d):
    ax.add_collection(LineCollection(_segmentos(G), colors="#d0d0d0", linewidths=.7))
    xs = [CIUDADES[n][1] for n in G.nodes]; ys = [CIUDADES[n][0] for n in G.nodes]
    ax.scatter(xs, ys, s=14, c="#bbbbbb", zorder=2)
    ex = res.orden
    ax.scatter([CIUDADES[n][1] for n in ex], [CIUDADES[n][0] for n in ex], s=38,
               c=color_expandidos, zorder=3, edgecolors="k", linewidths=.3)
    px, py = [CIUDADES[n][1] for n in res.camino], [CIUDADES[n][0] for n in res.camino]
    ax.plot(px, py, color="#d62728", lw=2.6, zorder=4)
    for n, c in ((o, "#2ca02c"), (d, "#ff7f0e")):
        ax.scatter(CIUDADES[n][1], CIUDADES[n][0], s=110, c=c, edgecolors="k", zorder=5)
        ax.annotate(n, (CIUDADES[n][1], CIUDADES[n][0]), xytext=(5, 5), textcoords="offset points",
                    fontsize=8, weight="bold")
    ax.set_title(f"{titulo}\nexpandidos: {res.expandidos}/{len(G)} nodos · coste {res.costo:,.0f} km",
                 fontsize=10)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect(1.1)


def dashboard_principal(G, filas_casos, datos, ruta_salida: Path):
    o, d = CASO_PRINCIPAL
    rd = next(f[2] for f in filas_casos if (f[0], f[1]) == CASO_PRINCIPAL)
    ra = next(f[3] for f in filas_casos if (f[0], f[1]) == CASO_PRINCIPAL)

    exp_d = np.array([x[4].expandidos for x in datos]); exp_a = np.array([x[5].expandidos for x in datos])
    t_d = np.array([x[4].tiempo_s for x in datos]); t_a = np.array([x[5].tiempo_s for x in datos])
    red = 100 * (1 - exp_a / exp_d)
    dist = np.array([x[2] for x in datos])

    fig = plt.figure(figsize=(18, 11.5))
    gs = GridSpec(3, 3, figure=fig, height_ratios=[1.35, 1, 0.9], hspace=.38, wspace=.22)

    dibujar_mapa(fig.add_subplot(gs[0, 0]), G, rd, "#1f77b4", f"DIJKSTRA  {o} → {d}", o, d)
    dibujar_mapa(fig.add_subplot(gs[0, 1]), G, ra, "#ff7f0e", f"A* (h = línea recta)  {o} → {d}", o, d)

    ax = fig.add_subplot(gs[0, 2]); ax.axis("off")
    kpis = [
        ("Reducción media de nodos expandidos", f"{red.mean():.1f} %", "#ff7f0e"),
        ("Reducción en el caso principal", f"{100*(1-ra.expandidos/rd.expandidos):.1f} %", "#ff7f0e"),
        ("Tiempo total Dijkstra / A* (>1 = A* más rápido)", f"{t_d.sum()/t_a.sum():.2f}x", "#1f77b4"),
        ("Pares con rutas distintas en coste", "0 de %d" % len(datos), "#2ca02c"),
    ]
    ax.set_title("INDICADORES CLAVE (1431 pares)", fontsize=11, weight="bold")
    for i, (txt, val, col) in enumerate(kpis):
        y = 0.86 - i * 0.24
        ax.text(0.02, y, val, fontsize=26, weight="bold", color=col, transform=ax.transAxes)
        ax.text(0.02, y - 0.08, txt, fontsize=9, color="#444", transform=ax.transAxes)

    ax = fig.add_subplot(gs[1, 0])
    etiquetas = [f"{a}\n→ {b}" for a, b, _, _ in filas_casos]
    x = np.arange(len(filas_casos)); w = .38
    ax.bar(x - w/2, [f[2].expandidos for f in filas_casos], w, label="Dijkstra", color="#1f77b4")
    ax.bar(x + w/2, [f[3].expandidos for f in filas_casos], w, label="A*", color="#ff7f0e")
    ax.set_xticks(x); ax.set_xticklabels(etiquetas, fontsize=7); ax.legend(fontsize=8)
    ax.set_title("Nodos expandidos por caso de estudio"); ax.grid(alpha=.3, axis="y")

    ax = fig.add_subplot(gs[1, 1])
    ax.hist(red, bins=20, color="#ff7f0e", edgecolor="k", alpha=.85)
    ax.axvline(red.mean(), color="k", ls="--", label=f"media {red.mean():.1f} %")
    ax.set_title("Reducción de nodos expandidos de A* vs Dijkstra\n(distribución en todos los pares)")
    ax.set_xlabel("% menos nodos expandidos"); ax.legend(fontsize=8); ax.grid(alpha=.3)

    ax = fig.add_subplot(gs[1, 2])
    ax.scatter(dist, exp_d, s=8, alpha=.35, c="#1f77b4", label="Dijkstra")
    ax.scatter(dist, exp_a, s=8, alpha=.35, c="#ff7f0e", label="A*")
    ax.set_xlabel("Distancia en línea recta origen–destino (km)"); ax.set_ylabel("Nodos expandidos")
    ax.set_title("Esfuerzo de búsqueda según la distancia"); ax.legend(fontsize=8); ax.grid(alpha=.3)

    ax = fig.add_subplot(gs[2, :]); ax.axis("off")
    cab = ["Origen → Destino", "Coste (km)", "Alg.", "Expandidos", "Generados", "Frontera máx.",
           "Tiempo (ms)", "Ahorro nodos"]
    celdas = []
    for a, b, r1, r2 in filas_casos:
        celdas.append([f"{a} → {b}", f"{r1.costo:,.1f}", "Dijkstra", r1.expandidos, r1.generados,
                       r1.max_abierta, f"{r1.tiempo_s*1000:.3f}", "—"])
        celdas.append(["", f"{r2.costo:,.1f}", "A*", r2.expandidos, r2.generados, r2.max_abierta,
                       f"{r2.tiempo_s*1000:.3f}", f"{100*(1-r2.expandidos/r1.expandidos):.1f} %"])
    tabla = ax.table(cellText=celdas, colLabels=cab, bbox=[0, 0, 1, 0.88], cellLoc="center")
    tabla.auto_set_font_size(False); tabla.set_fontsize(8)
    for (i, j), c in tabla.get_celld().items():
        if i == 0:
            c.set_facecolor("#2d1b69"); c.set_text_props(color="w", weight="bold")
        elif i % 2 == 0:
            c.set_facecolor("#fff1e0")
    ax.set_title("Tabla comparativa A* vs Dijkstra (mismo coste óptimo en todos los casos)",
                 fontsize=11, weight="bold", pad=4)

    fig.suptitle("DASHBOARD – Actividad 5: Ruteo GPS con A* y Dijkstra (red de %d ciudades de México)" % len(G),
                 fontsize=15, weight="bold", y=.995)
    fig.text(0.5, 0.005, "Equipo: " + " · ".join(EQUIPO), ha="center", fontsize=8, color="#555")
    fig.savefig(ruta_salida, dpi=140, bbox_inches="tight")


def dashboard_admisibilidad(G, razones_recta, razones_sin, razones_w, barrido, ruta_salida: Path):
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))

    meta = "Cancún"
    h_real = nx.single_source_dijkstra_path_length(G, meta, weight="weight")
    nodos = [n for n in G.nodes if n != meta]
    real = np.array([h_real[n] for n in nodos])
    ax[0].scatter(real, [h_recta(n, meta) for n in nodos], s=22, c="#1f77b4", label="h = línea recta")
    ax[0].scatter(real, [h_recta_sinuosa(n, meta) for n in nodos], s=22, c="#2ca02c", marker="s",
                  label="h = recta × 1.2")
    ax[0].scatter(real, [1.5 * h_recta_sinuosa(n, meta) for n in nodos], s=22, c="#d62728", marker="^",
                  label="h = 1.5 × recta × 1.2 (no admisible)")
    lim = real.max() * 1.05
    ax[0].plot([0, lim], [0, lim], "k--", lw=1, label="h = h* (frontera de admisibilidad)")
    ax[0].set_xlabel("Coste real mínimo h*(n) hasta Cancún (km)"); ax[0].set_ylabel("Heurística h(n) (km)")
    ax[0].set_title("Admisibilidad: h(n) ≤ h*(n) ⇔ por debajo de la diagonal"); ax[0].legend(fontsize=8)
    ax[0].grid(alpha=.3)

    ax[1].hist(razones_recta, bins=30, alpha=.7, color="#1f77b4", label="línea recta")
    ax[1].hist(razones_sin, bins=30, alpha=.7, color="#2ca02c", label="recta × 1.2")
    ax[1].hist(razones_w, bins=30, alpha=.5, color="#d62728", label="1.5 × recta × 1.2")
    ax[1].axvline(1, color="k", ls="--"); ax[1].set_xlabel("h(n) / h*(n)   (> 1 = sobreestima)")
    ax[1].set_title("Razón h/h* en todos los pares (nodo, destino)"); ax[1].legend(fontsize=8)
    ax[1].grid(alpha=.3)

    ws = sorted(barrido)
    exp = [np.mean(barrido[w]["exp"]) for w in ws]
    sub = [np.mean(barrido[w]["sub"]) for w in ws]
    ax[2].plot(ws, exp, "o-", color="#ff7f0e"); ax[2].set_xlabel("Peso w en f = g + w·h")
    ax[2].set_ylabel("Nodos expandidos (media)", color="#ff7f0e")
    ax2 = ax[2].twinx(); ax2.plot(ws, sub, "s--", color="#d62728")
    ax2.set_ylabel("Sobrecoste medio respecto al óptimo (%)", color="#d62728")
    ax[2].set_title("Weighted A*: menos nodos a cambio de perder optimalidad")
    ax[2].axvline(1, color="gray", ls=":"); ax[2].grid(alpha=.3)

    fig.suptitle("Análisis de admisibilidad de la heurística", fontsize=14, weight="bold")
    fig.tight_layout(); fig.savefig(ruta_salida, dpi=140, bbox_inches="tight")


# ----------------------------------------------------------------------------
# 8. PROGRAMA PRINCIPAL
# ----------------------------------------------------------------------------
def main():
    salida = Path(__file__).resolve().parent
    print("=" * 80)
    print("ACTIVIDAD 5 - Búsqueda con heurísticas: A* vs Dijkstra (ruteo GPS)")
    print("Equipo:", "; ".join(EQUIPO))
    print("=" * 80)

    G = construir_grafo(); adj = lista_adyacencia(G)
    print(f"Grafo: |V|={G.number_of_nodes()} ciudades, |E|={G.number_of_edges()} tramos "
          f"(k={K_VECINOS}, sinuosidad={FACTOR_SINUOSIDAD}) | conexo={nx.is_connected(G)}")
    print("Sucesores de CDMX:", [(v, round(w, 1)) for v, w in adj["CDMX"]], "\n")

    # --- Admisibilidad y consistencia (todas las metas posibles)
    print("VERIFICACIÓN DE ADMISIBILIDAD Y CONSISTENCIA (todos los pares nodo-destino)")
    candidatas = [("h = línea recta (Haversine)", h_recta),
                  ("h = recta x 1.2", h_recta_sinuosa),
                  ("h = 1.5 x recta x 1.2", h_ponderada(1.5)),
                  ("h = 3 x recta x 1.2", h_ponderada(3.0))]
    razones = {}
    print(f"{'Heurística':<30}{'Viol. admis.':>14}{'Viol. consist.':>16}{'max h/h*':>10}  Veredicto")
    for nombre, h in candidatas:
        va, vc, r, total = verificar_heuristica(G, h)
        razones[nombre] = r
        veredicto = "ADMISIBLE" if va == 0 else "NO admisible"
        print(f"{nombre:<30}{va:>10}/{total:<4}{vc:>16}{r.max():>10.3f}  {veredicto}")

    # --- Casos de estudio
    adjc = adj
    filas = estudio_casos(adjc)
    print("\nCASOS DE ESTUDIO (tiempo = mejor de 20 ejecuciones)")
    cab = f"{'Origen → Destino':<26}{'Alg.':<10}{'Coste km':>10}{'Expand.':>9}{'Gener.':>8}{'Frontera':>10}{'t (ms)':>9}{'Ahorro':>9}"
    print(cab); print("-" * len(cab))
    for o, d, r1, r2 in filas:
        for r in (r1, r2):
            ahorro = "" if r is r1 else f"{100*(1-r2.expandidos/r1.expandidos):.1f} %"
            print(f"{o + ' → ' + d:<26}{r.nombre:<10}{r.costo:>10,.1f}{r.expandidos:>9}"
                  f"{r.generados:>8}{r.max_abierta:>10}{r.tiempo_s*1000:>9.3f}{ahorro:>9}")
    o, d = CASO_PRINCIPAL
    rp = next(f for f in filas if (f[0], f[1]) == CASO_PRINCIPAL)
    print(f"\nRuta óptima {o} → {d} ({rp[3].costo:,.1f} km):\n  " + " → ".join(rp[3].camino))

    # --- Todos los pares
    print("\nEXPERIMENTO COMPLETO: todos los pares de ciudades (verificado contra NetworkX)")
    datos, barrido = todos_los_pares(G, adj)
    n = len(datos)
    e_d = np.array([x[4].expandidos for x in datos]); e_a = np.array([x[5].expandidos for x in datos])
    e_s = np.array([x[6].expandidos for x in datos])
    t_d = sum(x[4].tiempo_s for x in datos); t_a = sum(x[5].tiempo_s for x in datos)
    print(f"Pares evaluados: {n} | costo igual al óptimo (NetworkX) en Dijkstra y ambos A*: {n}/{n}")
    print(f"{'Algoritmo':<20}{'Expand. media':>15}{'Expand. mediana':>17}{'Generados medios':>18}")
    for nombre, idx in (("Dijkstra", 4), ("A* (recta)", 5), ("A* (recta·1.2)", 6)):
        ex = np.array([x[idx].expandidos for x in datos]); ge = np.array([x[idx].generados for x in datos])
        print(f"{nombre:<20}{ex.mean():>15.1f}{np.median(ex):>17.1f}{ge.mean():>18.1f}")
    print(f"Reducción media de nodos expandidos A*(recta) vs Dijkstra: {100*(1-e_a/e_d).mean():.1f} %")
    print(f"Reducción media A*(recta·1.2) vs Dijkstra               : {100*(1-e_s/e_d).mean():.1f} %")
    print(f"Tiempo total sobre los {n} pares: Dijkstra {t_d*1000:.1f} ms | A* {t_a*1000:.1f} ms "
          f"(aceleración {t_d/t_a:.2f}x)")

    print("\nWEIGHTED A* (f = g + w·h, h = recta x 1.2): efecto de perder la admisibilidad")
    print(f"{'w':>5}{'Expand. media':>15}{'Sobrecoste medio':>18}{'Sobrecoste máx.':>17}{'Pares óptimos':>15}")
    for w in sorted(barrido):
        b = barrido[w]
        print(f"{w:>5}{np.mean(b['exp']):>15.1f}{np.mean(b['sub']):>17.2f}%{np.max(b['sub']):>16.2f}%"
              f"{b['opt']:>10}/{n}")

    # --- CSV
    with open(salida / "a5_resultados.csv", "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["origen", "destino", "algoritmo", "costo_km", "expandidos", "generados",
                     "frontera_max", "tiempo_ms"])
        for o_, d_, r1, r2 in filas:
            for r in (r1, r2):
                wr.writerow([o_, d_, r.nombre, round(r.costo, 2), r.expandidos, r.generados,
                             r.max_abierta, round(r.tiempo_s * 1000, 4)])
        wr.writerow([]); wr.writerow(["weighted_A*", "w", "expand_media", "sobrecosto_medio_%", "sobrecosto_max_%"])
        for w in sorted(barrido):
            wr.writerow(["", w, round(np.mean(barrido[w]["exp"]), 2),
                         round(np.mean(barrido[w]["sub"]), 3), round(np.max(barrido[w]["sub"]), 3)])

    # --- Dashboards
    dashboard_principal(G, filas, datos, salida / "a5_dashboard.png")
    dashboard_admisibilidad(G, razones["h = línea recta (Haversine)"], razones["h = recta x 1.2"],
                            razones["h = 1.5 x recta x 1.2"], barrido, salida / "a5_admisibilidad.png")

    print("""
ANÁLISIS DE COMPLEJIDAD
  Dijkstra (heap binario): O((V + E) log V) tiempo, O(V) espacio. Explora en círculos
                           alrededor del origen sin saber dónde está la meta.
  A*                     : mismo peor caso O((V + E) log V), pero con una buena h expande
                           muchos menos nodos (se orienta hacia el destino). En el peor caso
                           (h = 0) es idéntico a Dijkstra. Espacio O(V).

CONCLUSIONES (contrástalas con los números impresos)
  - A* con h = línea recta devolvió el MISMO coste óptimo que Dijkstra en todos los pares,
    porque h es admisible (y consistente): la verificación numérica no halló violaciones.
  - A* expande menos nodos, pero en ESTE grafo (54 nodos) el tiempo de reloj es similar al
    de Dijkstra: calcular h cuesta y el grafo es muy pequeño. La ventaja en tiempo aparece
    en grafos grandes (miles o millones de nodos), donde expandir menos nodos pesa más.
  - El ahorro de nodos es mayor en rutas cortas y medias (CDMX → Monterrey: 58 %) y menor en
    rutas que cruzan todo el país (Tapachula → Nogales: 6 %), porque casi toda la red
    queda entre el origen y el destino y Dijkstra también debe explorarla.
  - Una heurística más informada (recta x 1.2, aún admisible) reduce más la búsqueda:
    cuanto más cerca de h* sin superarla, mejor.
  - Si h sobreestima (w > 1) se expanden menos nodos pero A* puede devolver rutas más
    largas: se pierde la garantía de optimalidad (ver Weighted A*).
  - Límite del modelo: la red es una aproximación y los pesos son distancias, no tiempos
    ni tráfico real.
""")
    print("Archivos: a5_dashboard.png, a5_admisibilidad.png, a5_resultados.csv")
    plt.show()


if __name__ == "__main__":
    main()
