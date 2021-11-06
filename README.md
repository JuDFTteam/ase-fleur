# ase-fleur
Package adding IO/calculator functionalities for the [FLEUR code](https://www.flapw.de) to the [ase](https://www.pypi.org/project/ase) package.

## IO functionality

By installing this package ase gains the ability to read/write the following files

- ``.xml`` files used by the Fleur code (read-only) (format name ``fleur-xml``)
- ``out.xml`` files produced by the Fleur code (read-only) (format name ``fleur-outxml``)
- input files for the Fleur input generator (``inpgen`` or ``input generator`` have to appear in the title comment for automatic detection) (format name ``fleur-inpgen``)

```python
from ase.io import read, write

atom = read('inp.xml')
atom = read('inp_fleur', format='fleur-inpgen')

write('new_inp', format='fleur-inpgen')
```

## Calculator

Setting up a Fleur executable:
```python
from ase_fleur import FleurProfile

profile = FleurProfile(['<path/to/fleur/executable>'],
                       ['<path/to/inpgen/executable>'])
```

Basic usage of Fleur calculator (is a very bare-bones implementation at the moment)

```python
from ase_fleur import Fleur

# Get an ASE structure to calculate
atoms = read(
    #...
    )

calc = Fleur(profile=profile)

calc.calculate(atom,
               properties=['energy'],
               [])
```
