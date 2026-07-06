document.addEventListener("DOMContentLoaded", function () {
  console.log("acre hero JS loaded");

  const hero = document.querySelector(".acre-hero");
  const img = document.getElementById("acre-hero-img");

  if (!hero || !img) {
    console.log("acre hero not found");
    return;
  }

  const svgWidth = 1512.00;

  const tileX0 = 10.5;
  const tileY0 = 54.64;

  const tileWidth = 300.21;
  const tileHeight = 332.59;

  const xStep = 300.21;
  const yStep = 430.90;

  const nCol = 5;
  const nRow = 10;

  /* 0.86 * 1.15 ≈ 0.99 */
  const zoom = 0.99;

  const positions = [];

  for (let row = 0; row < nRow; row++) {
    for (let col = 0; col < nCol; col++) {
      positions.push({
        x: tileX0 + col * xStep,
        y: tileY0 + row * yStep
      });
    }
  }

  let current = -1;

  function showPlot(index) {
    const pos = positions[index];

    const baseScale = hero.clientWidth / tileWidth;
    const scale = baseScale * zoom;

    const heroWidth = hero.clientWidth;
    const heroHeight = hero.clientHeight;

    const visibleTileWidth = tileWidth * scale;
    const visibleTileHeight = tileHeight * scale;

    const xPad = (heroWidth - visibleTileWidth) / 2;
    const yPad = (heroHeight - visibleTileHeight) / 2;

    img.style.width = `${svgWidth * scale}px`;
    img.style.transform = `translate(${xPad - pos.x * scale}px, ${yPad - pos.y * scale}px)`;

    current = index;
  }

  function randomPlot() {
    let next = current;

    while (next === current && positions.length > 1) {
      next = Math.floor(Math.random() * positions.length);
    }

    console.log("acre hero plot:", next);
    showPlot(next);
  }

  hero.addEventListener("click", randomPlot);

  window.addEventListener("resize", function () {
    if (current >= 0) showPlot(current);
  });

  randomPlot();
});

/* BIRD ANIMATION */
document.addEventListener("DOMContentLoaded", function () {
  console.log("acre birds JS loaded");

  if (window.matchMedia("(prefers-reduced-motion: reduce)").matches) {
    console.log("acre birds disabled: reduced motion");
    return;
  }

  if (document.querySelector(".acre-birds")) {
    console.log("acre birds already exist");
    return;
  }

  const stage = document.createElement("div");
  stage.className = "acre-birds";
  stage.setAttribute("aria-hidden", "true");
  document.body.appendChild(stage);

  for (let i = 0; i < 5; i++) {
    const perch = document.createElement("div");
    perch.className = "acre-bird-perch face-right";
    perch.dataset.birdIndex = String(i);

    const sprite = document.createElement("div");
    sprite.className = "acre-bird-sprite";

    perch.appendChild(sprite);
    stage.appendChild(perch);

    perch.addEventListener("mouseenter", function () {
      flyBird(perch, i);
    });

    perch.addEventListener("click", function () {
      flyBird(perch, i);
    });
  }

  console.log("acre birds created:", document.querySelectorAll(".acre-bird-perch").length);

  function rand(min, max) {
    return min + Math.random() * (max - min);
  }

  function clamp(x, min, max) {
    return Math.max(min, Math.min(max, x));
  }

  function setFacing(bird, dx) {
    /*
      Threshold avoids rapid left/right flickering when dx is tiny.
    */
    if (Math.abs(dx) < 0.7) return;

    if (dx < 0) {
      bird.classList.add("face-left");
      bird.classList.remove("face-right");
    } else {
      bird.classList.add("face-right");
      bird.classList.remove("face-left");
    }
  }

  function catmullRom(p0, p1, p2, p3, t) {
    /*
      Smooth curve interpolation through p1 -> p2.
    */
    const t2 = t * t;
    const t3 = t2 * t;

    return {
      x: 0.5 * (
        2 * p1.x +
        (-p0.x + p2.x) * t +
        (2 * p0.x - 5 * p1.x + 4 * p2.x - p3.x) * t2 +
        (-p0.x + 3 * p1.x - 3 * p2.x + p3.x) * t3
      ),
      y: 0.5 * (
        2 * p1.y +
        (-p0.y + p2.y) * t +
        (2 * p0.y - 5 * p1.y + 4 * p2.y - p3.y) * t2 +
        (-p0.y + 3 * p1.y - 3 * p2.y + p3.y) * t3
      ),
      s: 0.5 * (
        2 * p1.s +
        (-p0.s + p2.s) * t +
        (2 * p0.s - 5 * p1.s + 4 * p2.s - p3.s) * t2 +
        (-p0.s + 3 * p1.s - 3 * p2.s + p3.s) * t3
      )
    };
  }

  function samplePath(points, progress) {
    /*
      Duplicate endpoints so Catmull-Rom behaves nicely at start/end.
    */
    const p = [
      points[0],
      ...points,
      points[points.length - 1]
    ];

    const nSegments = points.length - 1;
    const scaled = clamp(progress, 0, 0.999999) * nSegments;
    const seg = Math.floor(scaled);
    const t = scaled - seg;

    return catmullRom(
      p[seg],
      p[seg + 1],
      p[seg + 2],
      p[seg + 3],
      t
    );
  }

  function smoothstep(t) {
    /*
      Gentle global acceleration/deceleration.
      Much less jerky than per-segment easing.
    */
    return t * t * (3 - 2 * t);
  }

  function buildFlightPath(bird, index) {
    const rect = bird.getBoundingClientRect();

    const startX = rect.left + rect.width / 2;
    const startY = rect.top + rect.height / 2;

    const vw = window.innerWidth;
    const vh = window.innerHeight;

    const minX = 70;
    const maxX = vw - 100;
    const minY = 70;
    const maxY = vh - 140;

    const p0 = { x: startX, y: startY, s: 1.00 };

    /*
      Give each bird a slightly different style:
      0 = broad sweep
      1 = high arc
      2 = low dip
      3 = lazy loop
      4 = fast crossing path
    */
    const style = Math.floor(rand(0, 5));

    const direction = Math.random() < 0.55 ? -1 : 1;
    const xBias = direction * rand(80, 180);

    /*
      First waypoint is now a real flight point, not a tiny vertical hop.
      This avoids the small stall immediately after takeoff.
    */
    let p1, p2, p3, p4, p5;

    p1 = {
      x: clamp(startX + direction * rand(140, 280), minX, maxX),
      y: clamp(startY - rand(120, 210), minY, maxY),
      s: rand(0.86, 0.96)
    };

    if (style === 0) {
      /* Broad smooth sweep across the page */
      p2 = {
        x: clamp(vw * rand(0.58, 0.80) + xBias, minX, maxX),
        y: clamp(vh * rand(0.35, 0.55), minY, maxY),
        s: rand(0.76, 0.88)
      };

      p3 = {
        x: clamp(vw * rand(0.28, 0.52), minX, maxX),
        y: clamp(vh * rand(0.14, 0.32), minY, maxY),
        s: rand(0.62, 0.76)
      };

      p4 = {
        x: clamp(vw * rand(0.12, 0.36), minX, maxX),
        y: clamp(vh * rand(0.42, 0.66), minY, maxY),
        s: rand(0.68, 0.80)
      };

      p5 = {
        x: clamp(vw * rand(0.42, 0.66), minX, maxX),
        y: clamp(vh * rand(0.26, 0.48), minY, maxY),
        s: rand(0.78, 0.90)
      };
    } else if (style === 1) {
      /* Higher, floatier arc */
      p2 = {
        x: clamp(vw * rand(0.50, 0.75), minX, maxX),
        y: clamp(vh * rand(0.16, 0.30), minY, maxY),
        s: rand(0.70, 0.84)
      };

      p3 = {
        x: clamp(vw * rand(0.20, 0.45), minX, maxX),
        y: clamp(vh * rand(0.08, 0.22), minY, maxY),
        s: rand(0.58, 0.72)
      };

      p4 = {
        x: clamp(vw * rand(0.34, 0.62), minX, maxX),
        y: clamp(vh * rand(0.28, 0.46), minY, maxY),
        s: rand(0.68, 0.82)
      };

      p5 = {
        x: clamp(startX - rand(160, 260), minX, maxX),
        y: clamp(startY - rand(145, 225), minY, maxY),
        s: rand(0.82, 0.94)
      };
    } else if (style === 2) {
      /* Dip downward, then recover */
      p2 = {
        x: clamp(vw * rand(0.58, 0.78), minX, maxX),
        y: clamp(vh * rand(0.48, 0.66), minY, maxY),
        s: rand(0.78, 0.90)
      };

      p3 = {
        x: clamp(vw * rand(0.36, 0.56), minX, maxX),
        y: clamp(vh * rand(0.62, 0.76), minY, maxY),
        s: rand(0.78, 0.90)
      };

      p4 = {
        x: clamp(vw * rand(0.18, 0.38), minX, maxX),
        y: clamp(vh * rand(0.26, 0.44), minY, maxY),
        s: rand(0.64, 0.78)
      };

      p5 = {
        x: clamp(vw * rand(0.44, 0.70), minX, maxX),
        y: clamp(vh * rand(0.24, 0.42), minY, maxY),
        s: rand(0.78, 0.90)
      };
    } else if (style === 3) {
      /* Lazy loop-ish path */
      p2 = {
        x: clamp(vw * rand(0.52, 0.70), minX, maxX),
        y: clamp(vh * rand(0.30, 0.48), minY, maxY),
        s: rand(0.74, 0.88)
      };

      p3 = {
        x: clamp(vw * rand(0.38, 0.58), minX, maxX),
        y: clamp(vh * rand(0.12, 0.26), minY, maxY),
        s: rand(0.60, 0.74)
      };

      p4 = {
        x: clamp(vw * rand(0.18, 0.34), minX, maxX),
        y: clamp(vh * rand(0.30, 0.52), minY, maxY),
        s: rand(0.66, 0.80)
      };

      p5 = {
        x: clamp(vw * rand(0.46, 0.68), minX, maxX),
        y: clamp(vh * rand(0.46, 0.62), minY, maxY),
        s: rand(0.78, 0.90)
      };
    } else {
      /* Faster crossing path */
      p2 = {
        x: clamp(vw * rand(0.68, 0.86), minX, maxX),
        y: clamp(vh * rand(0.26, 0.46), minY, maxY),
        s: rand(0.72, 0.84)
      };

      p3 = {
        x: clamp(vw * rand(0.26, 0.46), minX, maxX),
        y: clamp(vh * rand(0.20, 0.40), minY, maxY),
        s: rand(0.58, 0.72)
      };

      p4 = {
        x: clamp(vw * rand(0.12, 0.30), minX, maxX),
        y: clamp(vh * rand(0.44, 0.64), minY, maxY),
        s: rand(0.68, 0.82)
      };

      p5 = {
        x: clamp(vw * rand(0.50, 0.72), minX, maxX),
        y: clamp(vh * rand(0.34, 0.52), minY, maxY),
        s: rand(0.80, 0.92)
      };
    }

    /*
      Approach perch from above, with more random side offset.
    */
    const p6 = {
      x: clamp(startX - rand(70, 170), minX, maxX),
      y: clamp(startY - rand(110, 190), minY, maxY),
      s: rand(0.90, 0.98)
    };

    const p7 = { x: startX, y: startY, s: 1.00 };

    return [p0, p1, p2, p3, p4, p5, p6, p7].map(function (p) {
      return {
        x: p.x - startX,
        y: p.y - startY,
        s: p.s
      };
    });
  }

  function setFlapSpeed(bird, progress) {
    const sprite = bird.querySelector(".acre-bird-sprite");

    /*
      Takeoff and landing flap faster. Middle flight slower.
    */
    if (progress < 0.06) {
      bird.classList.add("is-taking-off");
      bird.classList.remove("is-landing");
      sprite.style.animationDuration = "0.30s";
    } else if (progress > 0.84) {
      bird.classList.remove("is-taking-off");
      bird.classList.add("is-landing");
      sprite.style.animationDuration = "0.34s";
    } else {
      bird.classList.remove("is-taking-off");
      bird.classList.remove("is-landing");
      sprite.style.animationDuration = "0.72s";
    }
  }

  function flyBird(bird, index) {
    if (bird.classList.contains("is-flying")) return;

    bird.classList.add("is-flying");
    bird.classList.add("is-taking-off");
    bird.classList.remove("is-landing");

    const points = buildFlightPath(bird, index);

    /*
      Much slower. This is the main knob.
      Increase more if still too quick.
    */
    const duration = rand(21000, 29000);

    let startTime = null;
    let lastPoint = points[0];
    let animationFrame = null;

    function cleanup() {
      if (animationFrame !== null) {
        cancelAnimationFrame(animationFrame);
      }

      bird.classList.remove("is-flying");
      bird.classList.remove("is-taking-off");
      bird.classList.remove("is-landing");
      bird.style.transform = "";

      const sprite = bird.querySelector(".acre-bird-sprite");
      sprite.style.animationDuration = "";
      sprite.style.backgroundPosition = "";

      setFacing(bird, 1);
    }

    function tick(timestamp) {
      if (startTime === null) startTime = timestamp;

      const rawProgress = clamp((timestamp - startTime) / duration, 0, 1);

      /*
        Global smoothing only. No per-turn braking.
      */
      const progress = smoothstep(rawProgress);

      const pos = samplePath(points, progress);

      const dx = pos.x - lastPoint.x;
      const dy = pos.y - lastPoint.y;

      setFacing(bird, dx);
      setFlapSpeed(bird, rawProgress);

      /*
        Small banking angle. Keep it subtle.
      */
      const bank = clamp(dx * 0.12 + dy * 0.015, -14, 14);

      bird.style.transform = `translate(${pos.x}px, ${pos.y}px) scale(${pos.s}) rotate(${bank}deg)`;

      lastPoint = pos;

      if (rawProgress < 1) {
        animationFrame = requestAnimationFrame(tick);
      } else {
        cleanup();
      }
    }

    animationFrame = requestAnimationFrame(tick);
  }
});
