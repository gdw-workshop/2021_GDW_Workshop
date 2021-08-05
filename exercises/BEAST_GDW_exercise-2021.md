Genomics of Disease in Wildlife BEAST Tutorial
==============================================


## First. All the files should be on your computer.  


The alignment we'll be analyzing today is the same as what we used in the Mr. Bayes tutorial yesterday.  Just as a reminder...the sequences are from two subtypes of a feline retrovirus (pumalentivirus A and B) that infect bobcats (Lru) and mountain lions (Pco).  We've constructed phylogenetic trees in BEAST using these and other sequence data from these viruses to understand host-pathogen relationships, intra-and inter-host transmission dynamics, and how the ecology of the viruses provide insight into the ecology of the hosts they infect.  Today we'll be doing a fairly standard BEAST workflow to estimate the phylogenetic relationships among these viral isolates.  The goal of this tutorial is to practice the steps required to perform an analysis using BEAST, and to introduce you to some of the concepts involved in the set up and interpetation of the analysis.


##  Second.  Some helpful resources for when you are trying this at home.  Some of these are more essential than helpful!


Keep the manual in hand and refer to it often when using BEAST!  It comes with the main BEAST package when you download the software.


There is a great BEAST user forum, which should be your first place to go when you run into a roadblock.  Chances are others have had the same issue so the answer may already be there but if not, you can post your question and someone will probably answer you soon.

https://groups.google.com/forum/#!forum/BEAST-users


There are great tutorials online for most of the analyses BEAST can perform so if you want to learn how to do more cool stuff with the program check them out: 

http://BEAST.bio.ed.ac.uk/tutorials

If you don’t want to bog down your computer for days to weeks when running BEAST on large datasets there is an open access server online.  All you need is your BEAST input file (which you are about to learn how to make!).

https://www.phylo.org/

Finally, there is a book from the creators of BEAST that provides background on many of the core functions of the software and goes into more detail about parameters, priors, etc.  I find it extremely useful to keep it by my side when using BEAST.

*Bayesian Evolutionary Analysis with BEAST* – Alexei J. Drummond and Remco R. Bouckaert (ISBN-13: 978-1107019652)



## Generating trees in BEAST is a multi-step process


Today we are going to (re)use an alignment of feline immunodeficiency virus isolates to build a dated tree in BEAST and along the way learn all of the additional programs that are necessary when using BEAST.

As a quick reminder from my overview earlier...

BEAUTi (Bayesian Evolutionary Analysis Utility) is a program that is used to create the input file for BEAST.  This is where you will select the type of analysis you want to run and the provide information about priors for the parameters that the analysis will be estimating.  The output of BEAUTi is an xml file that can be modified manually for customized BEAST runs.

BEAST (Bayesian Evolutionary Analysis Sampling Trees) uses a Bayesian MCMC algorithm to produce rooted phylogenetic trees and estimate many values for evolutionarily important parameters along the trees (node dates, evolutionary rates, etc.).

Tracer is a program used to evaluate .log files that are produced during a BEAST run so that you can evaluate the quality of the MCMC run and the quality of the output parameter estimates.  

TreeAnnotator will find the *best* tree from the thousands of trees sampled during the BEAST MCMC run and summarizes posterior parameter estimates on the tree.

FigTree is used to visualize the final tree and can be used to label nodes, branches, etc. for publication.



# Ok.  Lets get Started!

You should now have a *BEAST_Part1* folder in *GDW_Data* that has five files in it.  

Go to the *BEAST* folder in the *GDW-Apps* directory on your desktop.  All of the following programs are located here except for FigTree and Tracer which are in their own folders within GDW-Apps.


## Open BEAUTi v1.10.4

Once in BEAUTi...

**File > Import Data**

Now navigate to the directory where you put the files from Google Drive (Desktop/GDW_Data/BEAST_Part1)..
	Select *PLVAB_aln.nex*

You should see a summary of the file that you loaded under the *Partitions* tab.  Don’t worry about the details now.  If you don’t get an error message you’re doing well.

If you want to get complicated later on and analyze multiple data sets (i.e. different loci) or simultaneously analyze non-genetic traits you would import those additional data sets here before moving on.  Each data set would be imported as a different data partition so you would see multiple lines hear with each line describing a data partition to be analyzed.

The *Taxa* tab is primarily used when you want to estimate certain aspects of evolutionary history of pre-defined groups of related taxa.  We won’t use this today...but feel free to explore later :)

Select the *Tips* tab

Select *Use tip dates*

Click *Parse Dates*


Now click *Defined by a prefix and its order*

Under *Order* select *last*

Under *Prefix* select the *period* (.)

Select *Parse as a number*

Leave the rest of the settings at default (unchecked) and hit *OK*

The *Date* column of the *Tips* tab should now be populated with dates ranging from the 1980s to the 2010s.  Ignore the precision column for now. 

### Q1a: What did we just do?  Take a minute to think about it before moving on.  Why is this important in the types of analyses BEAST is typically used for?

TIP: Be sure to start the date with a different prefix if you do a full calendar date, which I now do when I can. For example Puma.1985-01-07. You can "Parse as a calendar date" as well but still good to start date with something else. 

### Q1b: What is the height column and why does that matter for analyses/output in BEAST?


Click on the *Traits* tab (even though we aren’t using it today:)).  

This is where you would define traits that you want to model over the tree.  
You could use location, phenotype, host association, serotype, vector species, etc.  
But this is beyond the scope of today’s introduction so let’s move on to the *Sites* tab.


### Ok.  This is where we start to get into the good stuff.


Here’s where you should ideally know enough about your sequences and the genes/organisms they came from that you can inform the model you’re using to analyze them.  
Today you will just have to take my word for it that the following parameters are good ones for these data.  When you are ready to analyze your own data in BEAST you can use information from the literature, model selection software like jModelTest, ore previous experience with your organisms of interest to select priors for your parameters.

Under *Substitution Model* select *HKY*

Use estimated base frequencies

Use a gamma distribution to estimate site heterogeneity (Q: What is this referring to?)

Its often good to keep the number of Gamma Categories low unless you have reason to do otherwise...select 4.

Partition codons into two groups *(1+2),3*

Use the first two check boxes to unlink the substitution rate and rate heterogeneity across the codon partitions we just defined.  This will allow different rates to be estimated for the (first/second) and (third) codon positions.  I tend not to estimate different base frequencies across the codon positions but the third box is up to you.

The last two boxes on this page - *Yang96* and *SRD06* - are options for 2 predefined sets of the above parameter options that you can use if you want.  We have basically used the *SRD06* set of parameters today.

Move on to the *Clocks* tab.

Under *Clock type* chose: *Uncorrelated relaxed clock*

The default *Relaxed Distribution* is Lognormal.  Leave this as is.


### Q2a: What do you think the difference is between a strict clock and relaxed clock?
### Q2b: What is an *uncorrelated* relaxed clock?
### Q2c: For what types of data sets (taxa, genes, evolutionary questions, etc) might one or the other of these options be most appropriate?



Next up: *Trees* tab!


For today we will keep it simple and use a *Constant Size Coalescent* Tree Prior but this is one area that BEAST has grown in recent versions.

These are coalescent priors, where the effective population size (Ne) varies through time according to a certain function (Ne(t)).  

The first set of models can reflect population size changes since the MRCA. This can be particularly important for pathogens and also reflects an estimation where you should have a decent idea of what type of growth your taxa has undergone. The Bayesian Skyline option can estimate the timeline of historical changes in effective population size based on the pattern of coalescence in your dataset. Look for these figures in the phylodynamics lecture this afternoon. 

Leave the default of *Random Starting Tree* and move on to the *States* tab.

If you want to do ancestral reconstruction you can specify some options here.  But we don’t today.  
So let’s move on to the *Priors* tab.

This is another place it is important (although not necessarily critical) that you have some knowledge of the sequences/taxa that you are analyzing.  Each of the lines on this page are present because of the parameter choices we have made thus far.  
If you select different model parameters next time you use BEAUTi...you will have different priors when you open this tab. 

## Please take a minute to read the description of each prior and seee if you can figure out how they relate to the model choices we've made in the *Sites*, *Clocks*, and *Trees* tabs. 

For today we will leave them all at default except the following:

*ucld.mean*: change this to lognormal with an intial value of 0.1

Click on the *Operators* tab.  Today (and generally) you don’t need to mess with these as long as the *Auto Optimize* box is checked in the upper left corner.  Sometimes the output from a run will give you a warning that how the chain samples a certain parameter needs to be tweaked and you can use this tab to do just that.  

Click on the *MCMC* tab.

Set the length of the chain to *1,000,000* and log the parameter estimates every *1,000*.  

These numbers are inadequate but it will let you get output files fast (and I have an output file from longer chains we can look at :)). 

[As a rule of thumb you want to end up with ~10,000 logged parameters/trees at the end of a BEAST run. This means that if you run the chain for 10^8 steps, you’ll log parameter estimates (and trees) every 10^4 steps. Make sense? If you run the chain less steps you log more frequently.  This is at least a good way to start.]

You can change the output file stem name if you want or leave it as is. If you start to do multiple combinations of BEAUTi parameters and BEAST runs from the same alignment file it is nice to label the output differently.  I use numbers, then letters, then obscenities in that order as I go through the process many times until it comes out right! 

Ok.  Click *Generate BEAST File* and save it where you want (the default is in the same location of your input alignment file).

A pop up window will appear.  This is your chance to change any additional parameters from default but today we will leave them alone.

Click Continue.

Wahoo!  Step one done!

NOTE: In your actual datasets you will want to use the Marginal likelihood estimation (MLE) after having run MCMC analyses that have converged on a posterior. These analyses can then start where that left off. 


## Open BEAST (also located in GDW-Apps/BEAST)


Load your newly created *file.xml* (PLVAB_aln.xml unless you changed it).  You can leave all of the check boxes and fields at default.

Select *Run* and you are off to the races.  

Easy compared to BEAUTi right?

It will take about 3 minutes to run this analysis.


[If you didn’t make it through BEAUTi or if you run into errors when you run BEAST, you already downloaded a BEAST input file that will work so go to your *GDW_Data/BEAST_Part1* directory and find *PLVAB_aln_GDW.xml* and use it to run BEAST.

Helpful tip...when the beast run completes, some important summary statistics that may be used to modidify subsequent runs are printed to the screen.  You can save the beast output as a .txt file to keep them for later reference.

While BEAST is running…


## Open Tracer (in Desktop/GDW-Apps/Tracer v1.7.1)


**File > Import Trace File** 

Select your *file_stem.log* file (in *GDW_Apps/BEAST_Part1*) that is in progress from your current BEAST run (you can view it before the BEAST run is complete).   

You should also open *PLVAB_aln_GDW1.log* which is a BEAST output file from the same alignment we used earlier but with a few extra parameter estimates and a chain that was run 100x as long.  It will work well to let you see what a log file will look like after the end of a sufficiently long MCMC chain.

You can also drag and drop files into the *Trace Files* area and open multiple trace files simultaneously to compare runs.

Look at the mean posterior estimate of each parameter value within the ‘Traces’ pane on the left...do they make sense?  

Ok. Probably none of this will make sense (because this is likely your first time doing this and you don’t know much about the dataset used) but this is where you will evaluate parameter estimates to make sure they make sense when you analyze your own data!  This is also the first set of results you've generated that can be informative for your study and can end up in your Nature publication.

This is also a good time to go back to the prior distributions and values we entered into BEAUTi to see how our choices may have influenced the outcome, or which parameter estimates may be wildy different than we thought they would be.

The effective sample size (*ESS*) value is an important metric to use to evaluate if you have enough samples from your chain to have accurate estimates of your parameters.  (Remember this from Mr. Bayes yesterday?) Before you do your own BEAST analyses...read up on ESS in one of the sources of knowledge I listed at the start of the tutorial.

If this is less than 100 it will be red to give you a warning that estimates should not be trusted.  

Between 100 and 200 they are yellow…caution.  

Above 200 is considered acceptable and there is likely little benefit to going beyond that.  

To increase low ESS values you can: 

1) run your chain longer (remember where you change this setting in BEAUTi?)

2) sample your chain more frequently

3) run multiple independent runs of BEAST and combine them (not covered here but you use *LogCombiner* to do this)

4) chose a less complex model as your data may not be informative enough for highly complex parameter estimates.

You can visualize the distribution/frequency of values sampled throughout the length of the MCMC chain by clicking on the tabs above the right pane.  

The *Trace* is especially useful to view in order to know if you have achieved good sampling of the posterior…it should look like a ‘spiny caterpillar’ (or at least that’s how it was taught to me).  

You can also select two parameters (by holding command) and then click on *Marginal Prob Distribution* to see how they relate to one another.  Try this for the CP1+2 and CP3 kappa values (Do you remember what kappa is? If not check back to the priors tab of Beauti or google it).  Do the same for CP1+2 and CP3 mu values.

## Q3: Were we correct to estimate values for these codon partitions separately?   What does this say about the flexibility of different codon positions to mutate/evolve over time?

Ok.  Spend as much time as you want with Tracer but when you’re ready let’s move on.


## Open TreeAnnotator (in Destktop/GDW-Apps/BEAST)



Go down toward the bottom and choose your *Input Tree File*, which is one of the files output from BEAST ending in *.trees*.  This file contains information about the trees that BEAST sampled throughout your MCMC run.
	
I have a file you can use here as well since your BEAST run may not be done yet and (it probably isn’t the best example anyway): *PLVAB_aln_GDW1.trees* . Be patient when going through these steps.  When I tested this out on the GDW laptops with the large input file it was taking a while to process things.


Choose your output file name ending with *.tre* and the location you want to save it. 


Now go back up and enter a value equal to 10% of the total number of MCMC steps in your chain (remember you entered this value in the *MCMC* tab in BEAUTi?) into the *Specify the burnin as the number of states* field.
	
i.e. if you are analyzing your BEAST output from a run of 10,000,000 steps then you’d have a burnin of 1,000,000.  For a chain length of 1,000,000 steps (like the inadequate example I had you do) enter 100,000.  
(Hint: If you forget the number of steps you have in a given chain you can open the log file in *Tracer* and it will tell you how many steps were in the run.)

Select *Maximum clade credibility tree* and *Median node heights*


Sweet.

Now select *Run* and let it go!  Don’t worry…It doesn’t take too long.


This program is annotating a single maximum clade credibility (MCC) tree (the tree with the highest node support values from all 10,000 or so trees logged in your .trees file) with lots of the posterior parameter values that were estimated across all of the 10,000 or so states in your .log file.  


Now let’s open the MCC tree and check it out.  This is the moment we’ve been working towards.  Exciting I know.


## Open FigTree (in Desktop/GDW_Apps/FigTree)


Use **File > Open** 
	Select your *.tre* file that you made in TreeAnnotator (or you can open PLVAB_aln_GDW1.tre)

Now you can use the panel on the left to play with how the tree looks, color branches and nodes, and visualize many of the parameters we estimated during the BEAST run.  

For example, click *Node Labels* and look through the options of what can be displayed then click on *Node Ages*.  Because we entered our sample dates in years (when making the input file in BEAUTi) the numbers now displayed on each node correspond to median posterior estimate of the age of that node (the number of years since the most recent sample that each node coalesced/diverged).  

You can also use the *Node Labels* to view the posterior support for each branch/node (same as Mr. Bayes...this value is similar to a bootstrap value on ML trees). 

Play around with the different options and see what you can learn about the samples and model estimates.

For extra credit you can open another window of FigTree with MrBayes output tree we constructed on Wed and compare the two.  While the branching patterns, nodes, and major clusters should be very similar between the two trees, you should notice a very major difference between them...hmmm.  What on Earth is going on???

Congratulations!  That’s your first run through BEAST.  If you’ve made it this far you deserve a pat on the back and a beer.


## Here are some brief answers to the questions you contemplated along the way...

A1a: BEAST produced time-trees when *Use Tip Dates* is selected.  A sample date (known or estimated) will incorporated into the evolutionary analysis at the same time the tree and other parameters are being evaluated by the MCMC chain.  The theory behind why this works is based on Kingman's Coalescent...the details of which are beyond the scope of this exercise but [here's](http://www.sfu.ca/biology/courses/bisc869/869_lectures/MHP_Coalescent.pdf) an easy to understand introduction to the concept if you're interested.

A1b: In most other types of trees, x-coordinate (branch length) of a sample is based on its genetic distance from the root.  In contrast, the x-axis of time-trees is basically a timeline starting on the left with current time (or the date of the most recent sample), and working backwards in time toward the root of the tree (the theoretical common ancestor to all of the sampled taxa).  Therefore, the *height* you see in the Tips tab is the difference in time units of each sample from the most recent sample.  Still confused?  I don't blame you.  Difficult to explain and difficult to understand.  Here's an example: if I'm analyzing a set of sequences collected between 2000 and 2018, the hight of my samples from 2018 will be zero and the hight of the samples from 2000 will be 18.  This will hopefully make more sense when we view the final tree later on.

A2a: A strict clock assumes all lineages on a tree (branches, taxa, etc) evolve at the same rate.  A relaxed molecular clock allows different lineages to have different rates. 

A2b: An uncorrelated relaxed clock means that even closely related lineages are allowed to have different rates because the rates between related branches/lineages/taxa are not correlated to one another.

A2c: This one is more of food for thought.  If you spend some time with the BEAST book and user forums described at the beginning of this tutorial you will start to gain an understanding of how this part all fits together :)

A3: When two parameters have posterior support for distinct distributions, this lends strong support that they behave differenty and therefore models that treat them differently are capturing more information about how these viruses evolve than models that treat all codon positions the same.  This is one example of looking at the parameter values of the posterior distribution to gain insight into the evolutionary process of your samples.


LAB BEAST part 2 – Confirming whether a molecular clock can be reliably estimated
=================================================================================

The shapes of phylogenetic trees reflect biological processes occurring on wide-ranging timescales, from as little as few weeks in rapidly evolving microbes to millions of years for macroevolutionary processes such as speciation. Regardless of the scale, being able to calibrate trees (i.e. to put a time-scale on them) and to date key events is an important part of understanding the dynamical processes that gave rise to the phylogenetic patterns we observe. Central to this idea of tree calibration and dating is the concept of a **molecular clock**: the stochastic but regular accumulation of nucleotide substitutions through time. Where these substitutions are accumulating fast enough to be picked up within sampling intervals, such populations are referred to as measurably evolving (Drummond et al. 2003) and molecular clocks can be directly estimated from genetic sequence data. 

Many disease-causing organisms, especially RNA viruses but also bacteria and some DNA viruses, have been shown to be measurably evolving (Biek et al, 2015). However, how do you know that this holds true for a given pathogen genomic data set?  **Before we are attempting to estimate a molecular clock from our data using BEAST, we should confirm that the data contain sufficient signal to do so.** The BEAST output in TRACER does not generally give us reliable indication in that respect: the program will always attempt to fit a molecular clock if we asked it to do so and will produce an estimate whether or not there is a genuine signal of measurable evolution in the data. 


## Testing for a molecular clock signal using TempEst

A standard way to check for ‘clock-like’ evolution in your data is to check whether there is an accumulation of substitutions over time, i.e. whether tree tips corresponding to samples collected later in time are more genetically divergent from the tree root (which represents the most recent common ancestor, mrca) than tips with earlier dates. What we need to test for such a pattern is a phylogenetic tree with dated tips. Importantly, this tree needs to have branches measured in genetic distance, not time, so it can’t be a tree estimated using a method already assuming a molecular clock (like BEAST). For this exercise, we will be looking at data set of genomes for the bacterial pathogen *Borrelia burgdorferi*, the agent of Lyme disease. 

Download the files for this exercise from *BEAST_Part2* from the github GDW/exercises/.  

You are provided with a Maximum Likelihood (ML) tree based on 23 *B. burgdorferi s.s.* genomes (main chromosome only, excluding data for the plasmids).  These samples were collected from a single site in Scotland during two time points, 1997 and 2013, so sixteen years apart. Because sampling was done at the same location both times, we can expect these isolates to be highly related. The key question is whether the isolates sampled at the later date are more genetically divergent from the root compared to the earlier isolates. In other words, **whether we see in an increase in root-to-tip divergence over time, consistent with a molecular clock**. 

The program we’ll use to test for a clock-like in the tree is called TempEst (Rambaut et al. 2016) and has been written by some of the same people who developed BEAST. Type 'tempest' in the terminal to start the program.  pen the program.  A java window will open...use it to navigate to *GDW_Data/BEAST_Part2* and import the ML tree of the *Borrelia* genomes (“Bbss_SubTree.nwk”). In the ‘Sample Dates’ tab you will see the names of the 23 isolates containing the year of sampling at the end. Similar to setting up a BEAST analysis in BEAUti, we need to extract the date information from the names using the ‘Guess Dates’ option. The dates are the last part of the name and the prefix is an underscore so specify this before hitting OK. Make sure the dates look correct in the Date column. Click the ‘Tree’ tab to take a look at the tree and then the ‘Root-to-tip’ tab to examine the pattern. 

*What would your conclusion be at this point – do you see evidence for measurable evolution and a molecular clock?*

It is worth considering what the root of the tree is that we are using as the reference point for calculating the root-to-tip divergence. TempEst expects a rooted tree and will use whatever root position it is provided with in the imported tree. In the current case, the tree had been previously rooted using an outgroup (other B burgdorferi s.s. genomes that fell outside of our clade of 23 genomes and that were subsequently removed). However, it is possible that we got the estimation of the root position wrong and that alternative root nodes should be considered. Clicking the ‘best-fitting root’ option on the top left provides a way to do this based on a range of criteria such as maximising the R2 of the correlation (all options produce the same outcome for our data but this is not always the case). 

*What would you conclude from the current plot? Can you give a first approximation of the molecular clock rate and the date of the most recent common ancestor?*

It is worth emphasising that the purpose of TempEst is not to estimate the molecular clock rate – you would still do this properly in program BEAST. But it gives a useful indication whether there is any clock-like signal in your sequence data, how much variation there is, and which sequences may be behaving very differently from the rest. If you don’t see a positive relationship between sampling date and genetic divergence at this stage, trying to estimate a clock rate in BEAST is probably not valid. 


## Examining the reliability of the molecular clock estimate obtained in BEAST

To save time, you won’t be asked to set up a BEAST analysis of the *Borrelia* data set and are instead provided with the log files of a finished run. Import the log file (“Bbss_chr_ucln.const.log”) into TRACER and examine the results for the clock rate. The analysis was done using a relaxed molecular clock model (Drummond et al. 2007)  and the parameter to look for is called ‘meanRate’. 

*How fast is B. burgdorferi evolving according to these results?*


As mentioned previously, BEAST will always attempt to estimate a clock rate. So how confident can we be that this rate estimate is really based on information contained in our data and not a spurious result? The most common way to assess this, is to **test whether similar rate estimates could have been obtained from the same data if the date information was not informative**. Randomising the dates of the sequences should achieve exactly that and is therefore a standard way to examine the reliability of the estimated clock rate. If we repeat the analysis with multiple data sets in which dates have been randomised and obtain values that are not overlapping with the original rate estimate, we can be confident that the originally estimated rate is a genuine product of our data and that it is statistically supported.

You are provided with the results from five such replicates with randomised dates (created with an R package called *TipDatingBeast*, which makes generating the alternative input files quick and easy). These log files are called ‘Bbss_chr_ucln.const.rep1.log’, ‘…rep2.log’, etc. Add them to the original log file in Tracer and examine the mean rate parameter. By highlighting all six files at the same time you can compare their highest posterior densities, which reflects the level of uncertainty around the mean parameter estimate. I would recommend using the ‘Marginal Density’ tab on the right for a visual comparison. 

*What do you conclude from these results regarding the clock rate of *B. burgdorferi*? Does the original rate estimate look distinct from its counterparts with randomised dates?*

The methods outlined here are generally applicable to any pathogen genomic data set (or subgenomic data). While most RNA viruses can be expected to be measurable evolving, especially if full genomes can be used, this is not always the case. Processes such as recombination for example can impact the relationship between sampling time and root-to-tip divergence in spurious ways. For double-stranded DNA viruses, Firth et al (2010) demonstrated that only some of the available genomic data sets exhibited measurable evolution. 

**Some general rules and recommendations:**

•	You should always assess the data for clock-like behaviour in TempEst before setting up a BEAST analysis. Even if the majority of the data exhibit the expected positive correlation, visual inspection is a good way to identify outliers. Such outliers might be attributable to errors in the data (such as a wrong date) or other processes (e.g. recombination). These outliers should be dealt with, either by correcting the data (where it is possible) or by removing them (if this can be justified). 

•	If the data look well behaved in TempEst, it should be fine to move ahead with a molecular clock analysis in BEAST. If they are not, or if the expected positive correlation can only be achieved under certain assumptions (e.g. best fitting root options), the BEAST analysis should be approached with caution. Confirming that the clock rate obtained is robust using the randomised date procedure or related techniques (Rieux and Khatchikian 2017) might be a good idea.

•	If the data don’t appear to contain sufficient temporal signal, several options could be considered. Expanding the time span of sampling for example (i.e. maximising the number of years between the oldest and most recent sample) might allow to reveal an increase in divergence over time, even if the data are generally noisy or the pathogen is evolving very slowly. Alternatively, more reliable clock rates might be available from previous studies of the same or related organisms and could be included as prior information in BEAST to guide the temporal calibration. 



**References**

Biek, R., Pybus, O. G., Lloyd-Smith, J. O., & Didelot, X. (2015). Measurably evolving pathogens in the genomic era. Trends in Ecology & Evolution, 30(6), 306–313. http://doi.org/10.1016/j.tree.2015.03.009

Drummond, A. J., Pybus, O. G., & Rambaut, A., Forsberg, R and A G Rodrigo (2003). Measurably evolving populations. Trends in Ecology & Evolution Vol.18 No.9, 481-488. http://doi.org/10.1016/S0169-5347(03)00216-7

Drummond, A. J., Ho, S. Y. W., Phillips, M. J., & Rambaut, A. (2006). Relaxed phylogenetics and dating with confidence. PLoS Biology, 4(5), 699–710. http://doi.org/10.1371/journal.pbio.0040088

Firth, C., Kitchen, A., Shapiro, B., Suchard, M. A., Holmes, E., & Rambaut, A. (2010). Using Time-Structured Data to Estimate Evolutionary Rates of Double-Stranded DNA Viruses. Molecular Biology and Evolution, 27(9), 2038–2051. http://doi.org/10.1093/molbev/msq088

Rambaut, Lam, de Carvalho & Pybus (2016). Exploring the temporal structure of heterochronous sequences using TempEst. Virus Evolution 2: vew007. http://dx.doi.org/10.1093/ve/vew007

Rieux, Adrien, and Camilo E. Khatchikian. (2017). "TipDatingBeast: An R package to assist the implementation of phylogenetic tip‐dating tests using BEAST." Molecular Ecology Resources 17, no. 4: 608-613. http://doi.org/10.1111/1755-0998.12603
