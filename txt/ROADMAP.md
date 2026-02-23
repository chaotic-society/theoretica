# Theoretica Development Roadmap

> **Last Updated:** 2026-02-23
> **Author:** M. Isgrò
> 
> This document outlines the strategic direction and major initiatives for the Theoretica project.

## Vision
The Theoretica project is an attempt to make powerful numerical methods available to the scientific computing community with an elegant and simple math-like interface. From battle-tested methods to cutting-edge research algorithms, Theoretica makes it possible to easily access numerical methods on a wide range of platforms. The proposed changes push Theoretica towards support for high-performance hardware, modern C++ standards, and integration with machine learning frameworks, making it a versatile and powerful tool for researchers and developers.

---

## 2026 Major Initiatives

### 🚀 Priority 1: Dual C++20 Adoption
**Status:** 🟡 **In Progress**
**Target:** Q1 2026  
**RFC:** [RFC002-CPP20-ADOPTION](https://github.com/chaotic-society/documents/blob/main/specification/theoretica/rfc/RFC-CPP20-ADOPTION.md)

**Brief Description:**
Considering the use cases of the library and its current strengths, Theoretica should strive to support both **high-performance hardware** and **embedded systems**. At the same time, Theoretica needs to modernize and adopt the C++20 standard, which is becoming the industry standard in scientific computing, to stay current and competitive.

Since embedded system support cannot always depend on C++20 compiler availability, a **hybrid, dual adoption strategy** is proposed:

- **Core modules** (linear algebra, real and complex analysis, basic utilities) will maintain **C++14 compatibility** with optional C++20 enhancements when available
- **Advanced modules** (autodiff, optimization, machine learning integration) will adopt **C++20** features fully

The `autodiff` module will particularly benefit from this change, as C++20 concepts provide clearer type constraints than C++14 SFINAE, aligning better with the library's design principles.


**Key Milestones:**
- [x] Define C++ standard separation strategy between core and advanced modules (Target: January 2026)
- [ ] Upgrade advanced modules to C++20 (`autodiff`, `optimization`) (Target: February 2026)
- [ ] Ensure and test C++14 compatibility of core modules (`algebra`, `complex`, `core`, `calculus`) (Target: February 2026)
- [ ] Document compiler requirements per module/header (Target: March 2026)
- [ ] Add embedded systems-specific features and testing (Target: March 2026)

**Technical Considerations:**
- Use preprocessor directives to detect C++ standard version
- Provide C++20 concept-based APIs alongside C++14 SFINAE when possible
- Maintain single codebase with conditional compilation for dual support

**Dependencies:**
- CI/CD updates for testing both C++14 and C++20 paths
- Updated build system (Make & CMake) configuration

---

### 🚀 Priority 2: Key missing features

**Status:** 🟡 **In Progress**
**Target:** Q1 2026
**RFC:** [RFC001-IO-MODULE](https://github.com/chaotic-society/documents/blob/main/specification/theoretica/rfc/RFC001-IO-MODULE.md)

**Brief Description:**
Some key features are currently missing from Theoretica that limit its usability and adoption. This initiative focuses on identifying and implementing these critical features, which may include:

1. **Comprehensive Test Coverage**
   - Complete test units for all modules using the Chebyshev framework
   - Ensure reliability and correctness across the entire library
   - Target: 90%+ code coverage for core modules
   - Target: 75%+ code coverage for advanced modules

2. **IO Module**
   - Facilitate data import/export in various formats
   - **CSV support**: Native implementation for lightweight, portable IO
   - **HDF5 support**: Integration for large-scale scientific data (optional dependency)
   - Create new `io` module with functions for reading and writing data structures to and from file

**Key Milestones:**
- [ ] Write missing test units for untested modules (Target: February 2026)
- [x] Implement CSV IO functionality in new `io` module (Target: February 2026)
- [ ] Implement HDF5 IO support (with optional build flag) (Target: March 2026)
- [ ] Document IO module with examples (Target: March 2026)

**Technical Considerations:**
- CSV: Native C++ implementation, no external dependencies
- HDF5: Optional dependency, enabled via CMake flag (`THEORETICA_ENABLE_HDF5`)
- API design should follow Theoretica's functional programming principles
- Integration with existing data structures such as `mat<>` and `vec<>`

**Dependencies:**
- **Chebyshev framework**: Feature-complete with context-based runtime and multi-threading support
- **HDF5 library**: External dependency for HDF5 support
- **Test data sets**: Sample CSV/HDF5 files for testing

---

### 🚀 Priority 3: GPU and CPU acceleration
**Status:** 🔴 Planning
**Target:** Q2 2026
**RFC:** Under development

**Brief Description:**  
Modern scientific computing software heavily relies on hardware acceleration to achieve high performance. Theoretica currently lacks acceleration features (except for the `parallel` namespace), which limits its adoption in high-performance computing scenarios. This initiative aims to research and implement GPU and CPU acceleration for key algorithms in Theoretica, leveraging modern libraries such as CUDA and SYCL for GPU support, OpenMP and possibly SIMD instructions for CPU optimizations. This will significantly enhance the performance of Theoretica in computationally intensive tasks, making it a more attractive option for researchers and engineers.

Technologies which should be considered for this initiative include:

| Technology | Maturity | Portability | Target Use Case |
|------------|----------|-------------|-----------------|
| **CUDA** | Mature | NVIDIA-only | GPU acceleration for NVIDIA hardware |
| **SYCL** | Modern | Cross-vendor (NVIDIA, AMD, Intel) | Portable GPU acceleration |
| **OpenMP** | Mature | Cross-platform | Multi-core CPU parallelism |
| **Highway** | Modern | Cross-platform | SIMD vectorization for CPU |

**Target Modules for Acceleration:**

1. **Linear Algebra** (`algebra` module)
    - Matrix operations
    - Matrix decompositions
    - Large-scale vector operations

2. **Numerical Integration** (`calculus` module)
    - Monte Carlo integration

3. **Signal Processing** (`signal` module)
    - FFT operations
    - Convolution

4. **Statistical Operations** (`statistics` module)
    - Large-scale sampling
    - Histogram generation

**Key Milestones:**
- [ ] Research and document GPU/CPU acceleration methods (Target: March 2026)
- [ ] Implement OpenMP support for multi-core CPU parallelism (Target: April 2026)
- [ ] Prototype CUDA acceleration for matrix operations (Target: May 2026)
- [ ] Implement SYCL support for cross-vendor GPU acceleration (Target: June 2026)
- [ ] Benchmark and optimize accelerated implementations (Target: August 2026)
- [ ] Document acceleration APIs and performance characteristics (Target: September 2026)

**Technical Considerations:**
- Maintain backward compatibility with CPU-only builds
- Minimize code duplication between CPU and GPU paths
- Handle memory transfers between CPU and GPU efficiently
- Document performance trade-offs for different problem sizes

**Dependencies:**
- CUDA Toolkit
- SYCL implementation (DPC++, ComputeCpp or hipSYCL)
- OpenMP compiler support
- Highway library (for SIMD)
- GPU hardware for testing and benchmarking

---

### 🚀 Priority 4: Machine learning integrations
**Status:** 🔴 Planning
**Target:** Q3 2026  
**RFC:** Under development

**Brief Description:**
Critical trends in machine learning should be researched and integrated into Theoretica, particularly **Physics-Informed Neural Networks (PINNs)** and **surrogate modeling**. Theoretica could benefit from integrating with popular ML frameworks to provide seamless interoperability, allowing users to leverage Theoretica's numerical methods within their machine learning workflows.

**Focus Areas:**

1. **Physics-Informed Neural Networks (PINNs)**
   - Solving differential equations using neural networks
   - Integration of Theoretica's automatic differentiation capabilities
   - Use of Theoretica's ODE/PDE solvers for validation

2. **Surrogate Modeling**
   - Use ML models to approximate expensive numerical simulations
   - Integration with optimization algorithms
   - Uncertainty quantification

3. **ML Framework Integration**
   - Python bindings for Theoretica (via pybind11 or similar)
   - Interoperability with TensorFlow, PyTorch, JAX
   - Export Theoretica data structures to ML frameworks

**Key Milestones:**
- [ ] Research PINNs, surrogate modeling, and current best practices (Target: June-July 2026)
- [ ] Implement Python bindings for core Theoretica modules (Target: July-August 2026)
- [ ] Develop prototype PINN implementation using Theoretica's autodiff (Target: September-October 2026)
- [ ] Create examples demonstrating ML integration (Target: November 2026)
- [ ] Major integration work and community feedback (Target: Q4 2026 - Q1 2027)


**Technical Considerations:**
- Python bindings may require exposing Python-friendly APIs
- Consider using `pybind11` or `nanobind` for Python interoperability
- Evaluate whether ML features should be separate library (`theoretica-ml`) or integrated
- Ensure automatic differentiation module (`autodiff`) is compatible with backpropagation

**Dependencies:**
- Python bindings library (pybind11, nanobind, or Cython)
- ML framework APIs (TensorFlow, PyTorch, or JAX)
- Interaction with the wider ML and scientific computing community
- Forward-mode auto-differentiation support in `autodiff` module if needed

---

### 🚀 Priority 5: Hybrid quantum computing
**Status:** 🔴 Planning
**Target:** Q4 2026  
**RFC:** Under development

**Brief Description:**
Quantum computing is an emerging field that holds promise for solving certain classes of problems more efficiently than classical computers. This initiative aims to explore the potential of hybrid quantum-classical computing within Theoretica. By researching current quantum computing frameworks (e.g., Qiskit, Cirq) and their applicability to scientific computing tasks, we can identify opportunities for integration. The goal is to enable Theoretica to leverage quantum computing resources for specific algorithms, enhancing its capabilities in cutting-edge scientific research.

Quantum computing is still an emerging technology, and this initiative will focus on research and exploration rather than immediate implementation. The primary objective is to understand how Theoretica can benefit from quantum computing advancements and prepare for future integration. In particular, current libraries and frameworks may not be ready for direct integration, but may be explored for future development.

**Important Note:** This is a **research and exploration initiative**, not a commitment to production-ready features. Quantum computing is still emerging, and this work focuses on understanding how Theoretica can position itself for future quantum computing advancements.

**Research Focus:**

1. **Quantum Computing Frameworks**
   - Evaluate Qiskit (IBM), Cirq (Google), Q# (Microsoft), Amazon Braket
   - Understand current capabilities and limitations (NISQ era)
   - Identify integration patterns

2. **Hybrid Quantum-Classical Algorithms**
   - Variational Quantum Eigensolver (VQE) for quantum chemistry
   - Quantum Approximate Optimization Algorithm (QAOA)
   - Quantum machine learning applications

3. **Applicability to Scientific Computing**
   - Linear algebra (quantum linear solvers like HHL algorithm)
   - Optimization problems
   - Monte Carlo sampling and integration
   - Quantum chemistry simulations

**Key Milestones:**
- [ ] Research hybrid quantum computing state-of-the-art (Target: October 2026)
- [ ] Evaluate quantum framework APIs and integration patterns (Target: November 2026)
- [ ] Prototype simple hybrid quantum-classical algorithm (if feasible) (Target: December 2026)
- [ ] Document findings and potential roadmap for quantum integration (Target: Q1 2027)
- [ ] Decision point: Continue, defer, or abandon quantum support (Target: Q1 2027)

**Technical Considerations:**
- Quantum computing may require Python bindings (most quantum SDKs are Python-based)
- NISQ-era quantum computers have significant noise and error rates
- Integration would likely be experimental/research-grade, not production-ready
- May be too early for practical integration; focus is on positioning for the future

**Dependencies:**
- Access to quantum computing simulators and/or hardware
- Quantum framework APIs (Qiskit, Cirq, etc.)
- Collaboration with quantum computing researchers
- Python bindings (if not already developed for ML integration)

---

## Quarterly Breakdown

### Q1 2026 (Jan - Mar)

**Focus:** Modernization and Foundation
- [ ] Dual C++20 adoption implementation
- [ ] IO module with CSV and HDF5 support
- [ ] **Release:** Version 26.Q1 including C++20 support and IO module

### Q2 2026 (Apr - Jun)

**Focus:** Performance and Acceleration
- [ ] OpenMP support for CPU acceleration
- [ ] Prototype SIMD optimization with Highway
- [ ] CUDA support for GPU acceleration
- [ ] SYCL support for GPU acceleration
- [ ] **Release:** Version 26.Q2 with OpenMP and experimental GPU acceleration support

### Q3 2026 (Jul - Sep)

**Focus:** Advanced Acceleration and ML Research
- [ ] Complete and stabilize SYCL and CUDA implementations
- [ ] Begin exploratory projects on machine learning integration
- [ ] Research PINNs and surrogate modeling
- [ ] **Release:** Version 26.Q3 with full GPU support (CUDA + SYCL)

### Q4 2026 (Oct - Dec)

**Focus:** Machine Learning Integration and Quantum Exploration
- [ ] Exploratory research on hybrid quantum computing
- [ ] **Release:** Version 26.Q4 with machine learning integrations

---

## Ongoing Initiatives

These activities continue throughout all quarters:

### Maintenance & Stability
- Bug fixes
- Code quality improvements

### Performance Optimization
- Continuous benchmarking
- Algorithm refinements
- Memory efficiency improvements
- Profiling and optimization

### Documentation & Community
- Tutorial and example expansion
- API reference updates
- Contributor onboarding and mentoring

### Build System & Infrastructure
- CI/CD improvements
- Multi-platform testing (Linux, Windows, macOS)

---

## Future Exploration (Beyond 2026)

Ideas under consideration for 2027 and beyond:

- **Native Machine Learning Module**: Dedicated module for ML algorithms implemented in C++
- **Distributed Computing**: MPI integration for cluster computing
- **Cloud Platform Support**: Integration with cloud computing platforms
- **Automatic Parallelization**: Compiler-level optimizations for parallelism
- **Symbolic Computation**: Computer algebra system capabilities
- **Quantum Computing Production Features**: If 2026 exploration is successful

---

## How to Contribute

We welcome contributions aligned with this roadmap! Please:

1. **Check Issues:** Visit [GitHub Issues](https://github.com/chaotic-society/theoretica/issues) to find tasks related to each initiative
2. **Join Discussions:** Participate in [GitHub Discussions](https://github.com/chaotic-society/theoretica/discussions) for roadmap feedback
3. **Read Contributing Guide:** Review [CONTRIBUTING.md](https://github.com/chaotic-society/theoretica/blob/master/txt/CONTRIBUTING.md) for code standards
4. **Contact Maintainers:** Reach out for guidance on where to start

### Priority Areas for Contributors

**High Impact, Good First Issues:**
- Writing test cases for existing modules
- Documentation improvements and examples
- Cross-platform testing and bug reports
- CSV and HDF5 IO implementation

**Advanced Contributions:**
- GPU acceleration implementations
- Python bindings development
- Machine learning integration prototypes
- Performance optimization

### Proposing New Features

Before working on major features:
1. Open a GitHub Discussion to gauge interest
2. Open a GitHub Issue outlining your proposal
3. Get feedback from maintainers and community
4. Proceed with implementation once approved

---

## Feedback & Updates

This roadmap is a living document. We review and update it:
- **Monthly:** Progress updates on active initiatives
- **Quarterly:** Milestone assessment and priority adjustments
- **Annually:** Strategic direction review

**Have feedback?** Open a discussion on [GitHub](https://github.com/chaotic-society/theoretica/discussions) or an [issue](https://github.com/chaotic-society/theoretica/issues).

---

## Legend

- 🔴 **Planning** - Research and RFC phase
- 🟡 **In Progress** - Active development
- 🟢 **Complete** - Implemented and released
- ⏸️ **On Hold** - Paused pending other work
- ❌ **Cancelled** - No longer pursuing

## Success Metrics

### Quantitative Metrics (2026 Goals)

- **Performance:** 10x+ speedup on GPU-accelerated operations vs. CPU baseline
- **Code Quality:** 90%+ test coverage for core modules
- **Community:** 50%+ increase in GitHub stars and active contributors
- **Releases:** Four quarterly releases on schedule (26.Q1 through 26.Q4)
- **Documentation:** Complete API reference and tutorials for all major features

### Qualitative Metrics

- **User Satisfaction:** Positive feedback on GitHub and social media
- **Industry Recognition:** Citations in academic papers; conference presentations
- **Developer Experience:** Contributors report smooth onboarding and clear documentation
- **Technical Leadership:** Theoretica recognized as modern, performant C++ scientific library

---
