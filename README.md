# Osmordred : Unified RDkit new descriptors in c++

Osmordred was inspired by the Dragon, Padel, Mordred and other toolkits to generate empirical molecular features
Our goal focus only on 0D,1D  and 2D molecular descriptors fused with rdkit backend at cpp level to get very fast computation in parallel if needed.

## Remark on reproductibility:

I spent quite time to implement a descent Information Content descriptor version based on the first paper from 1984 where Basak describes in detail his method https://doi.org/10.1016/B978-0-08-030156-3.50138-7.
So our implementation of Information Content is not 100% identical to Basak, Padel and Mordred but it follows the core Basak logic within RDkit where aromaticity is "specific".
This was indeed during this period that I also implement the Triplet features from Basak team. 

## Future:
Current version is around 10k lines of codes in only one file. 
It will be great to better integrate and refactor python bindings.
Additionally a list of other descriptors were added to produce now 3586 individual features.

## Speed 
This is fully parallelized. Lapack was selected of it speed specially on the SVD decompostion of symmetrical squared matrix instead of Eigen3 solvers. 
The Lapack can produce very small fluctuation for almost zeros Eigen values and affect very slighlty few descriptors.   

## Installation:

### Method 1 : from scratch create a new environment python 3.11
```
./setup_env.sh

./build.sh

conda activate osmordred

pip install dist/osmordred-0.2.0-cp311-cp311-macosx_15_0_arm64.whl --force-reinstall 
```

### Method 2 : include into your current environement (for the moment python 3.11 for RDKit 2023.9.3)

---
Mac Arm64
```
rm -rf _skbuild src/osmordred.egg-info dist
conda install boost==1.82.0 eigen lapack ninja python-build rdkit==2023.9.3 blas='*=*openblas' -c conda-forge
conda run python -m build
pip install dist/osmordred-0.2.0-cp311-cp311-macosx_15_0_arm64.whl --force-reinstall 
pip show osmordred   # normally you can see the installation in your environement 
```


---
Linux
```
rm -rf _skbuild src/osmordred.egg-info dist
conda install boost==1.82.0 eigen lapack ninja python-build rdkit==2023.9.3 blas='*=*mkl' -c conda-forge 
conda run python -m build
pip install dist/osmordred-0.2.0-cp311-cp311-linux_x86_64.whl --force-reinstall
pip show osmordred   # normally you can see the installation in your environement 
```

#### example for testing:
```
cd test
pip install tqdm
python test.py
```

#### note that attempting to `import osmordred` from python running in the skbuild directory does not work, as it imports the `osmordred` subdirectory.

### Known issue:
Such a complex molecule causes lag due to intensive computations : 
```
c12c3c4c5c1c1c6c7c2c2c8c3c3c9c4c4c%10c5c5c1c1c6c6c%11c7c2c2c7c8c3c3c8c9c4c4c9c%10c5c5c1c1c6c6c%11c2c2c7c3c3c8c4c4c9c5c1c1c6c2c3c41
```
## License:

BSD-3-Clause
