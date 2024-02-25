# DSQRM: Distribution-on-scalar Single-index Quantile Regression Model for Handling Tumor Heterogeneity

# Introduction
## Abstract
**This paper develops a distribution-on-scalar single-index quantile regression modeling framework to investigate the relationship between cancer imaging responses and scalar covariates of interest while tackling tumor heterogeneity. Conventional association analysis methods typically assume that the imaging responses are well-aligned after some preprocessing steps. However, this assumption is often violated in practice due to imaging heterogeneity. Although some distribution-based approaches are developed to deal with this heterogeneity, major challenges have been posted due to the nonlinear subspace formed by the distributional responses, the unknown nonlinear association structure, and the lack of statistical inference. Our method can successfully address all the challenges. We establish both estimation and inference procedures for the unknown functions in our model. The asymptotic properties of both estimation and inference procedures are systematically investigated. The finite-sample performance of our proposed method is assessed by using both Monte Carlo simulations and a real data example on brain cancer images from TCIA-GBM collection.**

## Framework
<details>
  <summary> Workflow </summary>
  ![Framework](/workflow.png)
  
  Fig 1. The workflow of our proposed association analysis framework.
</details>


## Software Interface
![GUI](interface.png)\
Fig 2. The interface of the software. 



# Software Usage 
## <mark>1. Input </mark>
### Covariates of interest
> choose a file of format .csv or .mat, containing the table of covariates of interest, with corresponding variable names. (not including the intercept)
- **.mat file:** containing a variable named "xdesign", which is a table of size (n, p0), with the first column as the sample name, and the column names as the variable names;
  - **n**: sample size 
  - **p0**: the number of covariates
- **.csv file:** the first column is the sample name, and the rest columns are the variables treated as the covariates of interests. \
  **<mark>Note: the continuous variables need to be normalized before being loaded.**


### Images
> Click the button **Load Image Data** to choose the images to be analyzed, and the format of this input depends on the choice <mark>**"Input Type of Images".**

#### Input type of images
- **whole images (default):** 3D or 2D images of whole brain, with the corresponding masks for tumor regions. 
    - a folder containing the images & corresponding tumor masks;
    - name of image file: SampleName + "image" + ".nii" / ".nii.gz" 
    - name of mask file: SampleName + 'mask' + ".nii" / ".nii.gz"
      - **ntype**: the number of tumor subtypes, which is the length of the unique nonzero values in a mask. 
- **image pixels:** extracted pixels of the tumor region for each sample.
  - a matlab data file (.mat) containing the pixels extracted from the tumor region (named "tumor_pixels"), & the ratios of each tumor subtype (named "sub_tumor_ratios");
  - **tumor_pixels:**
    -  a cell array of length n, each cell contains the pixels extracted from the tumor region; 
    -  a matrix of size (n, N), if the number of pixels in the tumor region are the same;
  - **sub_tumor_ratios:** a matrix of size (n, ntypes), where each row should be sum up to 1. 

### Design Matrix
> The final design matrix is a combination of the **covariates of interest** & **the tumor ratios of subtypes**, which is a matrix of size (n, p), where **p = p0 + ntype - 1.**


## 2. Settings
> - **Number of Grids:** m, the number of grids for the measurement of log-quantile density transformation.
> - **Target Quantile:** $\tau$, a scalar within the range of (0,1), denoting the target quantile level of the quantile regression. 

<details> 
<summary> <mark>Optional Input</mark> </summary>

#### Initials
> Click the button **Initials** to choose a matlab data file (.mat) containing the initial values for **the functional coefficients** $\beta(s)$ and **the link function** $g(\cdot)$ and its **first derivative** $\dot{g}(\cdot)$. \
The variable names should be
- **beta0**: a matrix of size (p, m)
- **g0**: a matrix of size (n, m)
- **dg0**: a matrix of size (n, m)

#### Bandwidth
> numerical values within (0,1), controlling the smoothness
- h1
- h2
- h3
</details>


### Output Folder
> Click the button **Output Folder** to choose a folder where the extracted distributional representations & the results of our model to be saved. \
The results will be saved in this folder by the name **"results.mat"**. 



## <mark> 3. Extract Distributional Representation </mark>
> Click the button **Extract Distributional Representation**, to have a visualization:
> - the whole brain with tumor segmentation, 
> - the extracted density functions, 
> - the log-quantile density functions.
 
**The sliders can be adjusted to visualize the representations of other samples.**


<details> 
<summary> Procedures </summary>

#### Input Type of Images = "whole images"
  - load the image & mask files for each sample;
  - extract the tumor pixels & the ratios of each subtype; 
  - obtain the density & log-quantile density (LQD) according to the extracted pixels;
  - display the original image of whole brain, with the annotation of subtypes of tumors; 
  - display the extracted densities & LQDs. 

#### Input Type of Images = "image pixels"
  - load the tumor pixels & the ratios of each subytpe; 
  - obtain the density & log-quantile density (LQD) according to the extracted pixels;
  - display the extracted densities & LQDs. 

</details>


## <mark> 4. Fit the Model </mark>
### (1) Estimation Procedures
- **Run**: click the button to start the algorithm and find the estimators
  - the indicator light turns <span style="color:green"> **green**</span> while the program is running
- **Pause**: click the button to pause the process, and click again to resume the process
  - the indicator light turns <span style='color:yellow'> **yellow** </span> while the program is paused
  - note: it takes a while for the program to pause
- **Stop**: click the button to stop the program
  - the indicator light turns <span style="color:red"> **red** </span>

#### Dislay the estimators 
> click the button the display the fitted functional coefficients and the link function. 

### (2) Inference Procedures
#### Settings
- Significance level: default = 0.05
- Number of Bootstrap: default = 200

#### Simultaneous Confidence Bands (SCB)
> Click the button to start the bootstrap procedure for constructing SCB for both functional coefficients and link function\
> user can choose if to display the SCB together with the estimators by checking the box at the bottom of the right panel 

#### Hypothesis Test
> Click the button to start the hypothesis test procedure for one of the functional coefficients
> - idx: denotes to conduct the hypothesis test on which covariate
> - p-value: the corresponding p-value after the bootstrap procedures


# An Example of GBM Study
## Folder Structure
> **examples**: this folder contains both input files and output results of this example. 

### - Input
   1. **images & masks**: '/examples/TCGA_flair_single_slice'
   2. **covariates of interest:** 'examples/TCGA_GBM_covariates.csv'
   3. **initial values** for the functional coefficients & link function: 'examples/initials.mat'
### - Output
  - saved in the matlab file "/examples/results.mat"
   



## Visualization of the Result
![GUI](interface_example.png)\
Fig 3. An example of the analysis on GBM dataset using the software.

