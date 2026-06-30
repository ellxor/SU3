# SU(3) Trajectories

**Instructions:**
- Download repo using `git clone https://github.com/ellxor/su3 --depth=1`
- Use the [parameters.h](parameters.h) file to set the parameters for both the model and integration.
- Build the program using the provided [build.sh](build.sh) script on Linux or macOS.
  - Requirements: `build-essential` package on Linux, or CommandLineTools on macOS (run `xcode-select --install`)
- Run the code `./trajectory_simulation`
- This outputs a bunch of `.txt` files to the selected output directory.
- Use/modify the [example notebook](graph.ipynb) to generate plots from these logs

Feel free to contact [me](mailto:ec294@st-andrews.ac.uk) if you have any questions :)

**Instructions for Windows:**
- Install [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) using `wsl --install` from PowerShell.
- Follow setup instructions (use default distribution, probably Ubuntu)
- Install packages needed: `sudo apt install build-essential git`
- Create directory: e.g. `mkdir local` or change local for any folder name, match this in parameters.h
- Follow instruction above to build and run.
- A Linux partition should appear in explorer on the left after WSL install (I think), and one can access files using this.

