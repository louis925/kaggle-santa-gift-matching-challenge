### Kaggle - Santa Gift Matching Challenge  

# Monte Carlo and Min-Cost Max-Flow Approach

Solutions to the 2018 Kaggle optimization competition - [Santa Gift Matching Challenge](https://www.kaggle.com/c/santa-gift-matching)

Authors: Louis Yang and Sky Huang

Methods: Monte Carlo, Min-Cost Max-Flow (MCMF), and Mixed-Integer Programming (MIP).

### Requirement
Max-Flow and MIP approaches require [Google Ortools](https://developers.google.com/optimization/).

#### Installation on Microsoft Visual Studio 2017
1. Download and unzip [OR-Tools for C++ (Visual Studio 2017)](https://github.com/google/or-tools/releases/download/v6.4/or-tools_VisualStudio2017-64bit_v6.4.4495.zip)
2. Create an empty C++ Visual Studio solution for one of the folder.
3. Change configuration to "Release" and platform to "x64".
4. In Solution Explorer, right click on the project > Properties.
5. In VC++ Directories, ([or-tools] is the directory where you unpacked the or-tools archive.)
	1. add "[or-tools]\include" to "Include Directories".
	2. add "[or-tools]\lib" to "Library Directories".	
6. In Linker > Input, add "ortools.lib" to "Additional Dependencies".

### Descriptions
* input - Download and upzip "gift_goodkids_v2.csv", "child_wishlist_v2.csv", and "sample_submission_random_v2.csv" to here from Kaggle.
* max_flow_mc - Monte Carlo with Min-Cost Max-Flow approach
* max_flow_merge_mc - Improved Monte Carlo with Min-Cost Max-Flow approach which merges groups
* max_flow_relax_all - Solve the complete relaxed problem as Min-Cost Max-Flow problem
* max_flow_search
* mc_swap - Original Monte Carlo swapping approach
* mip_mc - Monte Carlo selection with MIP optimization
* results - Result will be saved at here
* score - Grading function