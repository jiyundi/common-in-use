# Quick Guide

**Motivation:** I saved all of the codes that I used frequently for some time but later forgot. Once I needed to remind myself of these, I always wasted my time googling and viewing online pages. 

Thus, I created this repository to let me quickly direct to my used codes (Python, Linux, IDL, etc.) and also hope to help other people viewing this now. 

Please note that the codes in this repository are usually from others' work or in the public domain. Use caution if you'd like to refer to these.

## `Python` plotting quick look
```
import matplotlib.pyplot as plt
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,5))
im = ax.imshow(data, aspect='auto', cmap='viridis')
plt.colorbar(im, ax=ax)
plt.show()
```
```
import matplotlib.pyplot as plt
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(4,2)) # (width, height)
plt.subplots_adjust(hspace=0.3, wspace=0.6) # (height, width)
ax[0].hist(data[mask], bins=30)
plt.show()
```

## VS Code remote connection to Arizona HPC
0. (One-time, optional) Generate SSH key on Windows:
   ```
   ssh-keygen -t rsa -b 4096
   ```
   and press `Enter` 3 times. Now upload the public keys to HPC
   ```
   scp C:\Users\14477\.ssh\id_rsa.pub jiyundi@filexfer.hpc.arizona.edu:~/id_rsa.pub
   ```
   and write in `authorized_keys`
   ```
   ssh jiyundi@filexfer.hpc.arizona.edu "mkdir -p ~/.ssh && chmod 700 ~/.ssh && cat ~/id_rsa.pub >> ~/.ssh/authorized_keys && chmod 600 ~/.ssh/authorized_keys && rm ~/id_rsa.pub"
   ```
   The following doesn't guarantte a password-free remote connection in Step 9, but it's a worth to try -- Modify VS Code SSH Configuration file: (1) In VS Code, press Ctrl + Shift + P, choose "Remote-SSH: Open SSH Configuration File...", Open `C:\Users\14477\.ssh\config` and add the following
   ```
   # 1. Jumping board node (the node where you have previously configured the unencrypted key)
   Host hpc-jump
       HostName filexfer.hpc.arizona.edu
       User jiyundi
       IdentityFile C:\Users\14477\.ssh\id_rsa
   
   # 2. Dynamic calculation node template (connecting specific calculation nodes)
   Host hpc-node
       HostName r6u24n2.puma.hpc.arizona.edu
       User jiyundi
       IdentityFile C:\Users\14477\.ssh\id_rsa
       ProxyJump hpc-jump
   ```
1. (One-time) Go to UA [ServiceNow](https://uarizona.service-now.com/sp?id=kb_article_view&sysparm_article=KB0011701&sys_kb_id=a83f1b551b5dda103578773bdc4bcbea&spa=1)
2. (One-time) Installation Instructions > Instructions > Click the red `vpn.arizona.edu` > Log in and download. 
3. Connect to UA HPC VPN by typing `vpn.hpc.arizona.edu` and authenticate.
4. Use your local terminal (PowerShell), log on to `jiyundi@hpc.arizona.edu`, and authenticate.
5. You should be on HPC Gatekeeper (netid@gatekeeper ~). Now type `shell` and continue.
6. You should be on Puma or similar clusters (like `(puma) [netid@wentletrap ~] $`). Now apply for nodes and CPU sources with
   ```
   interactive -a my_supervisor_netid -n 4 -t 1:00:00
   ```
   where `my_supervisor_netid` should be replaced, `4` is the number of nodes applied for, and `1:00:00` is a 1-hour time slot for you.
7. The wait time of the step above may be up to 60 - 120 seconds. Then the terminal will show your node name such as `r6u24n2`. For example
   ```
   [netid@r6u24n2 ~]$
   ```
   Now you can type `hostname` to copy the output to your clipboard, making it easier to log in to VS Code.
8. Open a remote connection in VS Code. (Installation of Remote SSH is required)
9. Choose "+ Add New SSH Host" and then use what you copied in Step 7 to log in to your node 
   ```
   ssh jiyundi@r6u24n2.puma.hpc.arizona.edu
   ```
11. When you see "Open Folder" in the Explorer on the left of VS Code, that means you are all good to go!

## `rsync` Commands (Windows)
Remove `--dry-run` if you are ready.
```
D:\cwrsync_6.4.7\bin\rsync.exe -av --delete --info=progress2 --dry-run /cygdrive/f/runs_nautilus/ /cygdrive/d/_RschArchives/RSCH3/kl_github/runs_nautilus/
```

## Extract best-fit and corner plots from HPC
1. **Best-fit plots** (*Powershell*)
```
for ($i=3; $i -le 141; $i++) { $s="{0:D3}" -f $i; cp ".\Slit_$s\best_fit_spec.png" "review_best-fit\Slit_${s}_best_fit_spec.png" -ErrorAction SilentlyContinue }
```
2. **Corner plots** (*Powershell*)
```
for ($i=3; $i -le 141; $i++) { $s="{0:D3}" -f $i; cp ".\Slit_$s\corner_all.png" "review_corner\Slit_${s}_corner.png" -ErrorAction SilentlyContinue }
```
3. **Check if `post.txt` files exist in folders** (*Powershell*)
```
Get-ChildItem -Path "runs_20260710\Slit_*" -Directory | Where-Object { -not (Test-Path (Join-Path $_.FullName "post.txt")) } | ForEach-Object { $_.Name -replace '^Slit_', '' }
```
4. **Delete specific files** (*Powershell*)
```
Get-ChildItem -Path "Slit_*\post.txt" -Directory | Remove-Item -Recurse -Force
```

## `git` Commands (ASTR513)
Check and generate an SSH key
```
ls -al ~/.ssh
ssh-keygen -t ed25519 -C "jiyundi@gmail.com"
cat ~/.ssh/id_ed25519.pub (and copy to GitHub > Profile > Settings > SSH and GPG keys > New SSH key)
ssh -T git@github.com
```
Git pull & push
```
git clone git@github.com:ua-2025q3-astr501-513/hw2-jiyundi.git
git add .
git commit -m "Major updates"
git pull git@github.com:ua-2025q3-astr501-513/hw2-jiyundi.git main
git push origin main
```

## Commonly-used `Linux` Commands
* `rm -r`: Delete non-empty directory
* `tar`: Zipping --- `tar -cvzf new_file_name.tar.gz folder/`
* `tar`: Unzipping --- `tar -xvzf file_to_be_unzipped.tar.gz`
* `touch`: Make a new txt file
* `nano`: Edit a text file. `Cmd`+`X` to finish and save as.
* `ps`:  Process monitor
* `ps`:  Process monitor in monitor
* `echo`: Print. E.g., `echo "Today is $(date)"`
* `printf '%03d\n' 4`: print `004`
* `for` loop:
  ```
  for i in {0..4}; do python main.py --slitID $i; done
  ```
* `for` loop to rename files:
  ```
  for i in {0..4}; do printf '%0004d ' $i; done`: print `0000 0001 0002 0003
  for i in {1..100}; do mv "$i.txt" "$(printf '%03d' $i).txt"; done
  ```
* 只想升级 Chrome
  ```
  sudo apt update
  sudo apt install --only-upgrade google-chrome-stable
  ```
  这里只有 Chrome 会升级，不会升级别的软件。

## `Python` Format and Print
| **Python 3**                 | Appearance                   |
|:---------------------------- |:---------------------------- |
| `"{:5.2f}".format(-1.230)`   | `-1.23` (5 digits)           |
| `"{:5.2f}".format(1.2300)`   | `_1.23` (1 space at front)   |
| `"{:4d}".format(42)`         | `__42`  (2 spaces at front)  |
| `"{:04d}".format(42)`        | `0042`                       |

## `Python` Read CSV tables
Go to [`readcsv8.py`](./readcsv8.py). Read 8-column csv and ignore the first 3 lines (`"infile"=open('filename')`).

## `Python` Plot
Go to [`plot.py`](./plot.py)

## `Astropy` FITS header
* See header in console --- Go to [`print_header.py`](./print_header.py)
* Save header as TXT file --- Go to [`save_header_to_txt.py`](./save_header_to_txt.py)

## `Astropy` Read 1D spectrum
Go to [`read_spec1d.py`](./read_spec1d.py)

## `Astropy` Cosmological Calculations
* $d_A$ $({\rm kpc}/{''})$ --- Go to [`scale_calc.py`](./scale_calc.py) Separation in transverse/projected proper kpc per arcsec at redshift $z$. 

## `IDL` Specpro
Go to [`specpro.pro`](./specpro.pro)

Go to [`specpro_modified.pro`](./specpro_modified.pro) (Changes: L#2335-2337 --> 2335-2343; L#4053 --> 4060-4067)

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


## 本地监修 [jiyundi.github.io](https://jiyundi.github.io)
### 1️⃣ 安装前提
安装 Ruby（Windows: [RubyInstaller](https://rubyinstaller.org/downloads/), 约 1 GB）
接着：
```
gem install bundler jekyll
cd /D D:\_Archived\Websites\jiyundi.github.io
```
在项目根目录下新建一个文本文件，命名为`Gemfile`（记得删掉`.txt`后缀），在开头加入以下两行：
```
gem "bigdecimal"
gem "logger"
```
保存并关闭 Gemfile 文件运行：
```
bundle install
gem install bigdecimal
gem install logger
```

### 2️⃣ 启动cmd并运行：
```
cd /D D:\_Archived\Websites\jiyundi.github.io
bundle exec jekyll serve
```
运行成功后，命令行会显示网站预览地址。
