# Quick Guide

**Motivation:** I saved all of the codes that I ever used frequently for some time but later forgot. Once I needed to remind myself of these, I always wasted my time in googling and viewing online pages. 

Thus, I created this repository to let me quickly direct to my used codes (Python, Linux, IDL, etc.) and also hope to help other people viewing this now. 

Please note the codes in this repository are usually from others' work or in the public domain. Use caution if you'd like to refer to these.

## Commonly-used `Linux` Commands
### `tar`: Zip & Unzip
* Zipping --- `tar -cvzf new_file_name.tar.gz folder/`
* Unzipping --- `tar -xvzf file_to_be_unzipped.tar.gz`

## `Python` Format and Print
| **Python 3**                 | Appearance                           |
|:---------------------------- |:------------------------------------ |
| `"{:5.2f}".format(d=1.23)`   | `_1.23` (5 digits, 2-decimal float)  |
| `"{:5.2f}".format(d=-1.23)`  | `-1.23` (5 digits, 2-decimal float)  |
| `"{:4d}".format(d=42)`       | `__42`  (two spaces at front)        |

## `Python` Read CSV tables
[Go to `readcsv8.py`](./readcsv8.py) Read 8-column csv and ignore the first 3 lines (`"infile"=open('filename')`).

## `Python` Plot
[Go to `plot.py`](./plot.py)

## `Astropy` FITS header
* See header in console --- [Go to `print_header.py`](./print_header.py)
* Save header as TXT file --- [Go to `save_header_to_txt.py`](./save_header_to_txt.py)

## `Astropy` Read 1D spectrum
[Go to `read_spec1d.py`](./read_spec1d.py)

## `Astropy` Cosmological Calculations
* $d_A$ $({\rm kpc}/{''})$ --- [Go to `scale_calc.py`](./scale_calc.py) Separation in transverse/projected proper kpc per arcsec at redshift $z$. 

## `IDL` Specpro
[Go to `specpro.pro`](./specpro.pro)

