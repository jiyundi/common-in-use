# Quick Guide

**Motivation:** I saved all of the codes that I ever used frequently for some time but later forgot. Once I needed to remind myself of these, I always wasted my time in googling and viewing online pages. 

Thus, I created this repository to let me quickly direct to my used codes (Python, Linux, IDL, etc.) and also hope to help other people viewing this now. 

Please note the codes in this repository are usually from others' work or in the public domain. Use caution if you'd like to refer to these.

## Commonly-used `Linux` Commands
### `tar`: Zip & Unzip
* Zipping --- `tar -cvzf new_file_name.tar.gz folder/`
* Unzipping --- `tar -xvzf file_to_be_unzipped.tar.gz`

## `Python` Format and Print
| **Python 3**                 | Appearance                   |
|:---------------------------- |:---------------------------- |
| `"{:5.2f}".format(-1.230)`   | `-1.23` (5 digits)           |
| `"{:5.2f}".format(1.2300)`   | `_1.23` (1 space at front)   |
| `"{:4d}".format(42)`         | `__42`  (2 spaces at front)  |
| `"{:04d}".format(42)`        | `0042`                       |

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

## `ffmpeg` Conversion
### Video AVI --> MOV
```
ffmpeg -i "D:\Pr_Project_Materials\PR50\*.avi" -c:v qtrle -pix_fmt argb "D:\Pr_Project_Materials\PR50\output_animation.mov"
```
### Chrome WEBP --> PNG
```
Get-ChildItem "*.webp" | ForEach-Object {
$outputFile = $_.DirectoryName + "\" + [System.IO.Path]::GetFileNameWithoutExtension($_.Name) + ".png" 
ffmpeg -i $_.FullName $outputFile
}
```

## Video's 日本語 transcribe
```
pip install openai-whisper
```
```
whisper your_video.mp4 --language Japanese --task transcribe
```
Note: Whisper 使用的预训练模型由 OpenAI 提供，有多种尺寸（如 Tiny|39 MB、Base、Small、Medium、Large|1.51G），模型越大，精度越高，但占用的存储空间也更大。默认情况下，如果你没有指定模型，Whisper 会下载 Large 模型（1.51G），可以指定一个更小的模型来减少下载时间和存储占用。如果你需要释放空间，可以删除：`del C:\Users\14477\.cache\whisper\model-large.pt`


