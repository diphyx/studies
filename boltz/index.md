# Boltz

Trimming Adapters and Filtering Low-Quality Reads with Conda

### **Introduction to Boltz**

[Boltz](https://github.com/jwohlwend/boltz.git) 
## Introduction

Boltz is the state-of-the-art open-source model to predict biomolecular structures containing combinations of proteins, RNA, DNA, and other molecules. It also supports modified residues, covalent ligands and glycans, as well as conditioning the prediction on specified interaction pockets or contacts.

All the code and weights are provided under MIT license, making them freely available for both academic and commercial uses. For more information about the model, see our [technical report](https://doi.org/10.1101/2024.11.19.624167). To discuss updates, tools and applications join our [Slack channel](https://join.slack.com/t/boltz-community/shared_invite/zt-2zj7e077b-D1R9S3JVOolhv_NaMELgjQ).

## Inference

You can run inference using Boltz with:

```
boltz predict input_path --use_msa_server
```
### Suggested System Requirements
For optimal performance, we recommend using an AWS instance from the `g6` category. These instances are equipped with NVIDIA GPUs and provide more than 32 GB of memory, ensuring efficient processing for Boltz's computational requirements.


| **Resource**       | **Recommended Specification** |
|---------------------|-------------------------------|
| **Memory (RAM)**    | 64 GB                        |
| **CPU Cores**       | 16                           |
| **CPU Type**        | Intel Xeon or AMD EPYC       |
| **GPU Device**      | NVIDIA A100 or V100          |
| **GPU Memory**      | 40 GB                        |
| **Storage**         | SSD with at least 1 TB       |
| **Network**         | High-speed internet (1 Gbps) |

Boltz currently accepts three input formats:

1. Fasta file, for most use cases

2. A comprehensive YAML schema, for more complex use cases

3. A directory containing files of the above formats, for batched processing

To see all available options: `boltz predict --help` and for more information on these input formats, see our [prediction instructions](docs/prediction.md).
#### hint
This environment variable ensures that matrix multiplication operations in PyTorch use medium precision, which can improve performance on certain hardware configurations without significantly affecting accuracy. It is particularly useful when running Boltz on GPUs with limited precision support.
```bash
export TORCH_FORCE_FLOAT32_MATMUL_PRECISION=medium
``` 

example 
```bash
export TORCH_FORCE_FLOAT32_MATMUL_PRECISION=medium
boltz predict /volume/boltz_inputs/file.yaml \
  --use_msa_server \
  --out_dir /volume/boltz_output/ \
  --cache /volume/cache \
  --recycling_steps 10 \
  --diffusion_samples 5
```

## Boltz-2: 

Boltz 2 introduces significant enhancements over its predecessor, Boltz 1, particularly in terms of performance, scalability, and integration. Key improvements include:

* Improved Parallel Efficiency: Boltz 2 demonstrates superior performance on modern multi-core and multi-node architectures. It scales effectively on distributed systems, making it more suitable for large-scale simulations.

* Enhanced Solver Accuracy: The updated numerical solvers in Boltz-2 provide increased accuracy, particularly for simulations involving high-gradient fields and non-equilibrium dynamics.

* Expanded Model Library: Boltz-2 includes a broader set of pre-built physical models and boundary conditions, enabling researchers to simulate a wider variety of physical phenomena with minimal setup.

* Accelerated GPU Support: Boltz-2 provides optimized GPU acceleration, significantly reducing computation times for large datasets compared to Boltz-1.

![Boltz-2 Performance Comparison](boltz2_performance.png)
[Boltz-2 Performance Comparison](https://bit.ly/boltz2-pdf)

This figure clearly illustrates the performance differences between Boltz 1 and Boltz 2. It shows a comparative benchmark where Boltz 2 achieves over 2x speedup on a 64-core cluster, highlighting its superior scalability.

These upgrades make Boltz-2 a powerful tool for researchers requiring high-performance, accurate, and scalable simulations across a variety of scientific domains.
