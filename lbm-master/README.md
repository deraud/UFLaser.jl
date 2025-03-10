# Python Lattice Boltzmann 
Low level interface of Lattice Boltzmann technique as described in the PhD thesis of Nils Thuerey. Heavily depends on Numba for optimization.  

### Example: dam break
<img src="https://raw.githubusercontent.com/Maarten-vd-Sande/lbm/master/examples/dambreak.gif" width="400" height="400" />

### Example: cylinder flow
<img src="https://raw.githubusercontent.com/Maarten-vd-Sande/lbm/master/examples/cylinder_flow.gif" width="400" height="400" />

### Installation
To run the examples install an environment with conda:


```console
user@comp:~$ conda env create -f requirements.yaml
user@comp:~$ conda activate lbm
(lbm) user@comp:~$ python examples/cylinder_flow.py
```
