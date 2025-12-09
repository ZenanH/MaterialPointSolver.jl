```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: MaterialPointSolver.jl
  text: Material Point Method (MPM)
  tagline: A High-performance Backend-agnostic Material Point Method Solver in Julia
  actions:
    - theme: brand
      text: View on GitHub 👀
      link: https://github.com/LandslideSIM/MaterialPointSolver.jl
  image:
    src: /logo.png
    alt: MaterialPointSolver.jl

features:
  - icon: 🚀
    title: High-performance
    details: Carefully crafted kernel functions and implementation methods, leveraging the advantages of the Julia language for very fast computation speed!

  - icon: 🖇️
    title: Backend-agnostic
    details: A set of code can run on different hardware acceleration backends without worrying about data transmission issues.

  - icon: ⚽️
    title: MPM Algorithms
    details: The explicit solver comes with various built-in MPM-related algorithms, ready to use out of the box.

  - icon: 😎
    title: Plugin
    details: Your own algorithm? No problem, supports prototype design at any level, still fast and simple!
---
```