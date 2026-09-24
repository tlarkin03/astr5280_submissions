# ASTR 5820 — Origin and Evolution of Planetary Systems

Course code repository, Fall 2026. Prof. Paul O. Hayne, University of Colorado Boulder.

**Start with [SETUP.md](SETUP.md).** Thirty minutes, once, before Week 2.

`astrotools` is a single toolbox built up over the semester: each problem set adds
functions to it, and later sets import what earlier ones wrote. The tests in
`tests/` are the specification for each coding exercise — read them before writing
any code.

```bash
git clone https://github.com/phayne/astr5820.git
cd astr5820
conda env create -f environment.yml
conda activate astr5820
pip install -e .
python check_setup.py
```
