(function () {
  "use strict";

  var SRC = "./assets/data/dose-response-results.json";
  var FONT = "IBM Plex Mono";
  var BG = "#F3F5F7";
  var WASH = "#E7F1F0";
  var INK = "#1B2430";
  var MUTED = "#5C6773";
  var TEAL = "#0F6F6A";
  var VB = { w: 320, h: 180, l: 52, r: 300, t: 18, b: 142 };

  function logNorm(x, xmin, xmax) {
    return (Math.log(x) - Math.log(xmin)) / (Math.log(xmax) - Math.log(xmin));
  }

  function xOf(x, xmin, xmax, scale) {
    var t = scale === "log" ? logNorm(x, xmin, xmax) : (x - xmin) / (xmax - xmin);
    if (t < 0) t = 0;
    if (t > 1) t = 1;
    return VB.l + t * (VB.r - VB.l);
  }

  function yOf(y, ymin, ymax) {
    return VB.b - ((y - ymin) / (ymax - ymin)) * (VB.b - VB.t);
  }

  function el(name, attrs, text) {
    var node = document.createElementNS("http://www.w3.org/2000/svg", name);
    Object.keys(attrs).forEach(function (key) {
      node.setAttribute(key, attrs[key]);
    });
    if (text != null) node.appendChild(document.createTextNode(text));
    return node;
  }

  function tickLabel(v) {
    if (v >= 10) return String(Math.round(v));
    if (v >= 1) return String(Number(v.toFixed(2)));
    return String(Number(v.toFixed(2)));
  }

  function polyline(map, points) {
    return points
      .map(function (p) {
        return xOf(p.x, map.xMin, map.xMax, map.xScale).toFixed(1) + "," + yOf(p.y, map.yMin, map.yMax).toFixed(1);
      })
      .join(" ");
  }

  function draw(host, map) {
    var svg = el("svg", {
      class: "mc-curve__svg",
      viewBox: "0 0 " + VB.w + " " + VB.h,
      role: "img",
      "aria-label": map.aria || map.xLabel
    });
    svg.appendChild(el("rect", { width: String(VB.w), height: String(VB.h), fill: BG }));

    if (map.band) {
      var y0 = yOf(map.band.y1, map.yMin, map.yMax);
      var y1 = yOf(map.band.y0, map.yMin, map.yMax);
      svg.appendChild(el("rect", {
        x: String(VB.l),
        y: String(y0),
        width: String(VB.r - VB.l),
        height: String(y1 - y0),
        fill: WASH
      }));
      svg.appendChild(el("text", {
        x: String(VB.r - 4),
        y: String(y0 - 3),
        fill: TEAL,
        "font-size": "7",
        "font-family": FONT,
        "text-anchor": "end"
      }, map.band.label));
    }

    svg.appendChild(el("line", { x1: String(VB.l), y1: String(VB.t), x2: String(VB.l), y2: String(VB.b), stroke: INK, "stroke-width": "1" }));
    svg.appendChild(el("line", { x1: String(VB.l), y1: String(VB.b), x2: String(VB.r), y2: String(VB.b), stroke: INK, "stroke-width": "1" }));

    (map.yTicks || [map.yMin, map.yMax]).forEach(function (tick) {
      var y = yOf(tick, map.yMin, map.yMax);
      svg.appendChild(el("text", {
        x: String(VB.l - 6),
        y: String(y + 3),
        fill: MUTED,
        "font-size": "8",
        "font-family": FONT,
        "text-anchor": "end"
      }, tick.toFixed(1)));
    });

    (map.xTicks || [map.xMin, map.xMax]).forEach(function (tick) {
      var x = xOf(tick, map.xMin, map.xMax, map.xScale);
      svg.appendChild(el("text", {
        x: String(x),
        y: String(VB.b + 14),
        fill: MUTED,
        "font-size": "8",
        "font-family": FONT,
        "text-anchor": "middle"
      }, tickLabel(tick)));
    });

    svg.appendChild(el("text", {
      x: "14",
      y: "88",
      fill: MUTED,
      "font-size": "9",
      "font-family": FONT,
      transform: "rotate(-90 14 88)"
    }, map.yLabel));
    svg.appendChild(el("text", {
      x: "176",
      y: "172",
      fill: MUTED,
      "font-size": "9",
      "font-family": FONT,
      "text-anchor": "middle"
    }, map.xLabel));

    (map.series || []).forEach(function (series) {
      var color = series.color || (series.style === "dashed" ? INK : TEAL);
      var attrs = {
        fill: "none",
        stroke: color,
        "stroke-width": String(series.width || 2),
        "stroke-linejoin": "round",
        "stroke-linecap": "round",
        points: polyline(map, series.points)
      };
      if (series.style === "dashed") attrs["stroke-dasharray"] = "5 3";
      svg.appendChild(el("polyline", attrs));
      series.points.forEach(function (p) {
        svg.appendChild(el("circle", {
          cx: xOf(p.x, map.xMin, map.xMax, map.xScale).toFixed(1),
          cy: yOf(p.y, map.yMin, map.yMax).toFixed(1),
          r: "2",
          fill: color
        }));
      });
    });

    (map.callouts || []).forEach(function (note) {
      svg.appendChild(el("text", {
        x: xOf(note.x, map.xMin, map.xMax, map.xScale).toFixed(1),
        y: yOf(note.y, map.yMin, map.yMax).toFixed(1),
        fill: MUTED,
        "font-size": "7",
        "font-family": FONT
      }, note.text));
    });

    host.replaceChildren(svg);
  }

  function captionFor(fig, data, key) {
    var cap = fig.querySelector("figcaption");
    if (!cap) return;
    cap.setAttribute("data-dose-loaded", key);
  }

  function boot(data) {
    document.querySelectorAll("[data-dose-map]").forEach(function (host) {
      var key = host.getAttribute("data-dose-map");
      var map = data.maps[key];
      if (!map) return;
      draw(host, map);
      var fig = host.closest("figure");
      if (fig) captionFor(fig, data, key);
    });
  }

  function start() {
    fetch(SRC)
      .then(function (res) {
        if (!res.ok) throw new Error("dose-response-results.json " + res.status);
        return res.json();
      })
      .then(boot)
      .catch(function () {
        /* static SVG fallback already in the mount */
      });
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", start);
  else start();
})();
