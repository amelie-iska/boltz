<h1 align="center">Boltz-1:

Democratizing Biomolecular Interaction Modeling
</h1>

![](docs/image0.jpg)

Boltz-1 is an open-source model, named after Ludwig Boltzmann, which predicts the 3D structure of proteins, RNA, DNA and small molecules; it handles modified residues, covalent ligands and glycans, as well as the ability to condition the generation on pocket residues. Future improvements and finetunes of Boltz-1 may allow us to predict Boltzmann distributions of biomolecular complexes, similar to the way AlphaFlow, DiG, and MDGen do. Additionally, retraining Boltz-1 with a flow-matching objective instead of diffusion may allow us to speed up inference significantly, improve performance, and add additional functionality such as the ability to interpolate conformations or inpainting. 

The future of molecular modeling, what comes after Boltz-1, will include the ability to perform motif scaffolding and binder design with Boltzmann distributions rather than with single, static, low energy conformations. This will allow us to tune things like affinity and off rate, catalysis, and allosteric modulation. Suppose for example that we have a target protein in some low energy conformation, trapped inside of some energy well. Imagine being able to design a second protein with a specified Boltzmann distribution, where the new protein has a conformation with a much higher energy well, allowing the second designed protein to climb over the first energy barrier, and subsequently pull the first, target protein into some new conformation upon binding. 

When the second designed protein binds to the first, it introduces an interaction energy that modifies the first protein's energy landscape, effectively coupling their Boltzmann distributions. This binding lowers the energy barrier to a previously inaccessible conformation of the first protein, stabilizing a new target state. The first protein's Boltzmann distribution shifts to favor this conformation, as the interaction creates a new global minimum in the combined energy landscape. By designing the second protein's distribution to tune the interaction, researchers can precisely manipulate the thermodynamics and kinetics of the system, enabling controlled conformational transitions for functions like catalysis, affinity modulation, or allosteric regulation. This is the very near future of molecular engineering. 

For more information about the model, see the author's [technical report](https://gcorso.github.io/assets/boltz1.pdf).

## Installation
Install boltz with PyPI (recommended):

```
pip install boltz
```

or directly from GitHub for daily updates:

```
git clone https://github.com/jwohlwend/boltz.git
cd boltz; pip install -e .
```
> Note: we recommend installing boltz in a fresh python environment

## Inference

You can run inference using Boltz-1 with:

```
boltz predict input_path
```

Boltz currently accepts three input formats:

1. Fasta file, for most use cases

2. A comprehensive YAML schema, for more complex use cases

3. A directory containing files of the above formats, for batched processing

To see all available options: `boltz predict --help` and for more informaton on these input formats, see our [prediction instructions](docs/prediction.md).

## Training

If you're interested in retraining the model, see our [training instructions](docs/training.md).

## Contributing

We welcome external contributions and are eager to engage with the community. Connect with us on our [Slack channel](https://join.slack.com/t/boltz-community/shared_invite/zt-2uexwkemv-Tqt9E747hVkE0VOWlgOcIw) to discuss advancements, share insights, and foster collaboration around Boltz-1.

## Coming very soon

- [x] Auto-generated MSAs using MMseqs2
- [x] More examples
- [ ] Support for custom paired MSA
- [ ] Confidence model checkpoint
- [ ] Pocket conditioning support
- [ ] Full data processing pipeline
- [ ] Colab notebook for inference
- [ ] Kernel integration

## License

Our model and code are released under MIT License, and can be freely used for both academic and commercial purposes.


## Cite

If you use this code or the models in your research, please cite the following paper:

```bibtex
@article{wohlwend2024boltz1,
  author = {Wohlwend, Jeremy and Corso, Gabriele and Passaro, Saro and Reveiz, Mateo and Leidal, Ken and Swiderski, Wojtek and Portnoi, Tally and Chinn, Itamar and Silterra, Jacob and Jaakkola, Tommi and Barzilay, Regina},
  title = {Boltz-1: Democratizing Biomolecular Interaction Modeling},
  year = {2024},
  doi = {10.1101/2024.11.19.624167},
  journal = {bioRxiv}
}
```

In addition if you use the automatic MSA generation, please cite:

```bibtex
@article{mirdita2022colabfold,
  title={ColabFold: making protein folding accessible to all},
  author={Mirdita, Milot and Sch{\"u}tze, Konstantin and Moriwaki, Yoshitaka and Heo, Lim and Ovchinnikov, Sergey and Steinegger, Martin},
  journal={Nature methods},
  year={2022},
}
```
