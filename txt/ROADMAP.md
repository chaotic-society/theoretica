# Theoretica Development Roadmap

> **Last Updated:** 2026-09-19
> 
> **Author:** M. Isgrò
> 
> This document outlines the strategic direction and major initiatives for the Theoretica project.

## Vision
Theoretica is an open source project with the objective of making powerful and cutting-edge algorithms more accessible and transparent to the scientific computing community. This is done by implementing a wide range of algorithms, keeping up with new numerical methods and techniques, as well as providing access to these with a simple interface and a readable and documented codebase. Current initiatives push Theoretica towards supporting hardware acceleration and high-performance computing, modern C++ standards and possible extensions supporting relevant applications (e.g. machine learning).

### Objectives
The current objective of the project is to finalize current features and ensure stability and quality of the library, with the objective of delivering a first official release. This first release may still be in beta, with possible API changes dictated by user feedback. After these changes, a stable release is expected.

### How can I participate?
Contributions are welcome and appreciated, if you intend to contribute in any way to the project please have a look at the [related paragraph](https://github.com/chaotic-society/theoretica/edit/master/txt/ROADMAP.md#Contributing) or at the [Contributing guide](https://github.com/chaotic-society/theoretica/blob/master/txt/CONTRIBUTING.md).

---

## Initiatives

### Priority 1: Test Coverage

**Status:** 🟡 In Progress (82%)

**Issue:** [#53](https://github.com/chaotic-society/theoretica/issues/53)

Reliable scientific computing needs extensive testing of all algorithms and implementations. For this reason, a major task in Theoretica is writing test units for all of the library's modules. The Chebyshev framework, purposefully created for this task, is currently being employed. Most modules of the library have a good amount of tests, and at the moment only the `pseudorandom` and `statistics` modules are lacking proper test units. This is an ongoing initiative, as strengthening existing tests is also important and valuable.

#### Milestones
- [ ] Write base test units for all modules (82%)
- [ ] Extend existing test units

#### Technical Considerations
Test units are written using Chebyshev, which is documented in [TESTING.md](https://github.com/chaotic-society/theoretica/blob/master/txt/TESTING.md).

#### Dependencies
There are ongoing API changes in the `pseudorandom` module which may stall testing of some features.

---

### Priority 2: IO Module

**Status:** 🟢 Complete

**Issue:** [#100](https://github.com/chaotic-society/theoretica/issues/100)

Implement a new IO module supporting common operations for standard output and file IO, including support for CSV files and the HDF5 format.
The HDF5 format is a structured file format which allows IO of large-scale scientific datasets, and will not be implemented standalone.
The CSV file format is simple and straightforward to implement, so a standalone implementation will be included.

#### Milestones
- [x] Implement utilities for standard output
- [x] Implement support for CSV file format
- [x] Integrate HDF5 external library
- [x] Implement support for HDF5 file format

#### Technical Considerations
The new external library needs appropriate integration in the Makefile and CMakefiles, to allow enabling or disabling it and proper linking.

#### Dependencies
- External HDF5 library
- Test data sets for CSV and HDF5 formats

---

### Priority 3: C++20 Dual Adoption

**Status:** 🟡 In Progress (50%)

**Issue:** [#101](https://github.com/chaotic-society/theoretica/issues/101)

As C++20 is becoming the industry standard for high-performance scientific computing, Theoretica should move in the direction of adopting its most useful features, particularly `concepts`. A major concern in this sense is that Theoretica could have relevant applications in embedded systems, and fully adopting C++20 could hinder support for these platforms. For this reason, a dual adoption plan is being followed, with select features moving or being implemented in C++20, while the rest of the library remains compatible with C++14, with possible enhancements which get automatically turned on when the newer standard is supported. This is a partial solution for the beta release, as user feedback is fundamental in determining if a full upgrade to C++20 is desirable.

The `autodiff` module could particularly benefit from concepts for simplifying SFINAE and template resolution.

#### Milestones
- [x] Upgrade key modules/features to C++20 (`autodiff`, `optimization`)
- [x] Test C++14 compatibility
- [ ] Document C++20 requirements for upgraded features
- [ ] Deliberate on embedded systems support

#### Technical Considerations
- Auto-detect C++ standard with preprocessor directives
- Implement C++20 concepts as alternatives to C++14 SFINAE
- Both C++ standards should correctly compile and pass tests

#### Dependencies
- Add C++ standard configuration in Makefile and CMakeFiles
- Extended test matrix for C++14 and C++20 builds

---

### Priority 4: Hardware Acceleration

**Status:** 🟡 In Progress

Modern scientific computing relies on hardware acceleration on CPUs and GPUs to achieve high performance. A fully-featured, advanced scientific computing library must support hardware acceleration for both consumer hardware and HPC clusters. Support for parallelization using OpenMP is being implemented in the `parallel` module, while possible approaches for SIMD vectorization are being considered (SIMD intrinsics vs wrappers vs external libraries). Support for GPU acceleration should also be taken into consideration, especially for specific solvers and algorithms which could greatly benefit from it. An additional technology which should be taken into account for future support is OpenMPI, fundamental for HPC clusters. This could mean providing wrappers, schedulers and solvers using this platform.

Technologies under considerations include:

| Technology | Maturity | Portability | Target Use Case |
|------------|----------|-------------|-----------------|
| **CUDA** | Mature | NVIDIA-only | GPU acceleration for NVIDIA hardware |
| **SYCL** | Modern | Cross-vendor (NVIDIA, AMD, Intel) | Portable GPU acceleration |
| **Kokkos** | Modern | Cross-vendor | Portable CPU and GPU acceleration |
| **OpenMP** | Mature | Cross-platform | Multi-core CPU parallelism |
| **Highway** | Modern | Cross-platform | SIMD vectorization for CPU |
| **OpenMPI** | Mature | Cross-platform | Distributed HPC |

Specific features which could greatly gain from hardware acceleration:
- Large-scale linear algebra operations
- Monte Carlo methods
- Fourier transform and convolution
- Large-scale distribution sampling

#### Milestones
- [ ] Research and document acceleration methods
- [ ] Implement OpenMP parallelization for key features in `parallel` module
- [ ] Explore GPU acceleration with different technologies
- [ ] Benchmark and optimize accelerated implementations
- [ ] Document acceleration methods and hardware/software requirements

#### Technical Considerations
- Maintain compatibility with sequential execution of all algorithms
- Minimize code duplication between backends (OpenMP, SIMD, GPU)
- Run benchmarks to document trade-offs for specific algorithms

#### Dependencies
- CUDA Toolkit
- SYCL (DPC++, ComputeCpp, hipSYCL)
- Highway
- Various CPU and GPU hardware

---

### Proposals
The following are proposals for extensions, add-ons and branching projects for Theoretica which could be considered for implementation.

#### Proposal 1: Sparse linear algebra and PDEs
A key task in scientific computing is to solve PDEs, which usually entails using sparse matrices, linear algebra solvers and pre-conditioners.
Support for sparse linear algebra as well as key algorithms could make Theoretica a powerful toolbox for PDE solving.
This would include both the building blocks for solving PDEs as well as targeted solvers, potentially shipped in a dedicated add-on module.

#### Proposal 2: Machine Learning Integrations
Research machine learning applications of the library and consider integration possibilities with other frameworks and libraries for the task.
In particular, Physics-Informed Neural Networks (PINNs) and surrogate models could be interesting fields of application.
Integrating with other frameworks could entail developing Python bindings, at least for key modules.

---

## Contributing

We welcome contributions aligned with this roadmap!

1. **Check Issues:** Visit [GitHub Issues](https://github.com/chaotic-society/theoretica/issues) to find tasks related to each initiative
2. **Join Discussions:** Participate in [GitHub Discussions](https://github.com/chaotic-society/theoretica/discussions) for roadmap feedback
3. **Read Contributing Guide:** Review [CONTRIBUTING.md](https://github.com/chaotic-society/theoretica/blob/master/txt/CONTRIBUTING.md) for a complete guide on contributing to this project
4. **Contact Maintainers:** Reach out for guidance on where to start and general support. We are here to help!

### Key Contributions

**High Impact, Good First Issues:**
- Writing test cases for existing modules
- Documentation improvements and examples
- Cross-platform testing and bug reports

**Advanced Contributions:**
- New numerical methods
- CPU/GPU accelerated implementations
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

