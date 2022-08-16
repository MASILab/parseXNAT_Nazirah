To run on ACCRE




**First time installation**
```
module load Anaconda3
conda create --name mrtrix3
conda install -n mrtrix3 -c mrtrix3 mrtrix3
```

**Add to ~/.bashrc**

```
module load GCC/8.2.0
module load OpenMPI/3.1.4
module load FSL
module load MATLAB
module load Anaconda3
PATH="/home/mohdkhn/.conda/envs/mrtrix3/bin:$PATH"
```

**Update mrtrix3**

```
conda update -n mrtrix3 -c mrtrix3 mrtrix3
```
