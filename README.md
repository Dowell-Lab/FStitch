# Fast Read Stitcher (FStitch)
Fast Stitch Reader (FStitch) rapidly processes read coverage files into contigs of actie and inactive regions of transcription with its intended being primarily for refining annotations in nascent transcription data (GRO-seq, PRO-seq, NET-seq, etc.)<sup>1</sup>. Using FStitch, you can:

* Discover unannotated large regions of active transcription for enahncer and other non-coding RNA discovery
* Better annotate 5' and 3' ends of genes for differential transcription analysis and Tfit RNAPII modeling
* Filter regions of active transcriptional activity for downstream application (i.e. Tfit bidirectional predictions for enhancer identification)
* Differentiate and analyze genome-wide coverage of active transcription between treatment types

The following is an example of FStitch output:

![Alt text](images/IGV_SNAP.png)

*Integrative Genomics Viewer (IGV) snap shot demonstrates the annotations obtained using FStitch. Color ‘green’ indicates regions of inactive transcription (signal is not singificantly above background "noise"). Color ‘blue’ represents active transcription on the forward (pos) strand and ‘red’ on the reverse (neg) strand.*

## General Usage
Here are the minimal commands needed to run FStitch from start to finish; for greater detail on usage and file types see below. 
```
$ FStitch train -b </path/to/sample.bedGraph> -s (+/-) -t </path/to/TrainingFile>  -o </path/to/Parameters.hmminfo>

$ FStitch segment -b </path/to/sample.bedGraph> -s (+/-) -p </path/to/Parameters.hmminfo> -o </path/to/segmentFile.bed>
```

## System Requirements
FStitch is written in the C++ programming language, with C++11 support and uses OpenMP<sup>4</sup> to support multi-threading.  With this in mind, users will need to have a GCC compilers later than version 4.7 and < 7.1.0 to compile and run FStitch. To check you compiler version, 
```
$ gcc —-version
$ g++ —-version
```
Note, for those running FStitch on a compute cluster, commonly you will need to perform a ‘module load gcc<version>’ to compile FStitch. Please ask your sys admins for questions on module load behavior.
    
### Setup

Download the FastReadStitcher/ directory from this url or clone to your local machine. If your compiler is up to date, you can compile FStitch by moving into the FastReadStitcher/ directory and running 
```
$ sh setup.sh
=========================================
Sucessfully Compiled
```
In short, the setup.sh just runs “make clean” and "make" in the src/ directory. If everything compiles, you should see "Sucessfully Compiled" at the end. Importantly, you will now see the executable “FStitch” in the src directory.

## Running FStitch

The Fast Read Stitcher program is divided into two main commands: `train` and `segment` 
* `train`: estimates the necessary probabilistic model parameters
* `segment`: uses the output model parameters and segments the entire genome into *active* and *inactive* regions of transcription

To view the usage statement, use following standard arguments (must be compiled!):
```
$/src/FStitch -h 

$/src/FStitch --help
```

### Coverage File (bedGraph) Format Requirements

The provided coverage file must be in bedGraph format (See the UCSC description<sup>2</sup>) which is a BED4 file where the fourth column represents coverage over the annotated start/end positions. For example:

```
chr    start    end    coverage
1      0        100     3
1      107      117     1
```

***IMPORTANT***

Your .bedGraph file **should not contain 0 values** and should be **non-normalized**. FStitch performs an internal normalization.

There are two main tools for generating bedGraph coverage files from BAM files, deepTools and BEDTools. By default, deepTools bamCoverage will have discrete bins (50bp) and will therefore calculate average coverage over regions, rather that contigs of regions with equal read coverage, and "smooth" the data. While this is not a problem for visualizaiton at smaller bins, it will conflict with normalization. Therefore, we recommend using default BEDtools<sup>3</sup> genomecov settings:
    
```
$ bedtools genomecov -ibam <file.bam> -g <file.bedGraph> -bg -s <+/->
```

Specifying five prime (-5 argument) in the “genomecov” may allow for cleaner annotations however unspecified five prime bed works just fine as well. 

Often times, it is useful to merge the positive and negative strands for conversion to TDF if you are using IGV as your genome browser so you can view both strands on one track. While FStitch will only train and segment one strand at a time, you can provide this concatenated pos/neg bedGraph file as long as it is formatted with appropriate integers in the fourth column (i.e. -\<integer> for a negative strand coverage region). FStitch will remove regions based on which ever stand you specify in the -s/--strand argument. An example of how to generate an appropriately formatted concatenated pos/neg bedGraph file is as follows:

```
awk 'BEGIN{FS=OFS="\t"} {$4=-$4}1' ROOTNAME.neg.bedGraph \
 > ROOTNAME.neg.formatted.bedGraph

cat \
 ROOTNAME.pos.bedGraph \
 <(grep -v '^@' ROOTNAME.neg.formatted.bedGraph) \
 | sortBed \
 > ROOTNAME.sort.cat.bedGraph
```

Note that the last command, sortBed, will require BEDTools. This is not required for FStitch, but is good practice as it will process quicker and is required for conversion to TDF if so desired.


## FStitch train
FStitch uses two probabilistic models to classify regions of high read density into "active" transcriptional regions: Logistic Regression (LR) and a Hidden Markov Model (HMM). The LR coefficients are estimated via a user defined label training file. Because we are classifying regions as signal or noise, FStitch requires regions of the genome that show characteristic transcription or high read dense profiles and regions of the genome that display noise or not a profile of nascent transcription or a read dense region. With this information, FStitch trains a logistic regression classifier and then couples it to a HMM. The transition parameters for the Markov model are learned via the Baum Welch algorithm and thus do not require user labeled training data.  

In short, FStitch requires regions the user considers active transcription and regions considered inactive (noise). The more regions provided to FStitch the more accurate the classifications however we have in Cross Validation<sup>1</sup> analysis that roughly 15-20 regions of active and inactive regions will yield accurate classifications. 

### Training file Format

These regions are provided to FStitch using a specific file format with four columns separated by tabs: chromosome, genomic coordinate start, genomic coordinate stop, (0 if “noise” or 1 “signal”). An example is given below:

```
~/TrainingFile.bed
chr    start    end     label
1      1000     6000    1 
1      6001     8250    0
13     13500    19000   0
13     21000    23000   1
```

The segments do not need to be in any order and can be from any chromosome, however ***each region must not overlap any other segment*** as this will cause confusion in the learning algorithms for the LR classifier. This makes sense -- calling a region as both "ON" and "OFF" in binary world is confusing. 

**Very important**: If FStitch is being used on stranded data, the bedGraph file used in the `train` must correspond to the strand indicated in the *TrainingFile.bed*. For example, if the strand in the training file comes from the forward strand but the user supplies a bedGraph file that is on the reverse strand, then learned parameters will not be accurate. 

Running FStitch train is simple once you have your coverage data in the correct format and have created the training file above. The following is a description of arguments:

**Required Arguments**

|Flag|Type|Desription|
|----|----|----------|
|-b  --bedgraph| \</path/to/BedGraphFile>              | bedGraph File from above
|-s  --strand  | \<+/->                                | Specifes which strand (pos/neg) you trained on <br> **NOTE: You can only train on ONE strand!**</br>
|-t  --train   | \</path/to/TrainingFile.bed>          | Training File from above (BED4 format)
|-o  --output  | \</path/to/outDir/Parameters.hmminfo> | Training Parameter OutFile (.hmminfo extension)

**Optional Arguments**

|Flag|Type|Desription|
|----|----|----------|
|-n  --threads | \<integer>                            | number of processors, default 1

An example command is therefore:

    $ FStitch train -b </path/to/BedGraphFile.bedGraph> -s (+/-) -t </path/to/TrainingFile.bed> -o </path/to/Parameters.hmminfo>

The output will be a *Parameters.hmminfo* file that will store the learned parameters for the LR and the HMM transition parameters needed in `segment`. Here's an example of what the training file *Parameters.hmminfo* should look like.

```
####################################################
#                  Fast Read Stitcher
#Parameter Estimation Output
#Command Line                    :~/FastReadStitcher/src/FStitch train -s + -b bedgraphs/SRR000001.bedGraph -t fstitchTrain.bed -o Dowell2018.hmminfo
#Date/Time                       :2018-10-26 12:13:16
#Learning Rate                   :0.400000
#Max Iterations                  :100.000000
#Convergence Threshold           :0.001000
#####################################################
Converged                        :True
Final LogLikelihood              :-3530.215862
Logistic Regression Coefficients :1.013783,0.000000,-0.001297,186.575328
HMM Transition Parameter         :0.933049,0.066951,0.021232,0.978768
~
```

### Training Tips

Annotating regions of "active/singal" and "inactive/noise" transcription ("ON" and "OFF") and generating a desirable training file typically requires a little bit of trial-and-error. That said, the following are some tips to expedite the learning process.

##### Avoid "overtraining"
If your aim is to segment your sample into long contigs (e.g. for improving 5' and 3' gene annotation), it is best **not to overtrain**, or in other words provide too many training examples. If you provide >40 regions, you will get noticably larger .bed output files as a result of increasing the number of segments, or "contigs". In other words, the more regions you prodide, the more you simply obtain regions where coverage = 0 is "OFF" and any coverage > 0 is "ON". We don't need FStitch to give us this information.

##### Try to pick large training regions
If you annotate training regions that are <1000bp, there will likely not be enough annotated regions of coverage in the given bedGraph file to accurately calculate the LR coefficients. Try to pick a mixture of regions between 1b and 50kb for training.

##### Do not pick "dead zones"
For regions you annotate as "OFF", be sure to include regions that have some background signal, or "noise". Remember, you will be training using your bedGraph file, so if you pick regions that have "0" coverage and you followed the instructions above for coverage file formatting, you will not be providing FStitch any regions to assess. FStitch will try to account for data "gapiness" or the distance between read coverage regions, and therefore the more "background" (think of ChIP input controls) you can provide it, the better FStich will be able to distinguish the "noise" from the "activity".

## FStitch segment
FStitch `segment` uses the parameters obtained from `train` (from above, \</path/to/Parameters.hmminfo>) as input, as well as the original bedGraph file. A description of the arguments are given below.

**Required Arguments**

|Flag|Type|Desription|
|----|----|----------|
|-b  --bedgraph     | \</path/to/sample.bedGraph>           |BedGraph File Format from above
|-s  --strand       | \<+/->                                |Specifes which strand (pos/neg) you wish to segment <br> **NOTE: You can only segment on ONE strand at a time!**</br>
|-p  --params       | \</path/to/Parameters.hmminfo>        |Training parameters prduced from train module
|-o  --output       | \</path/to/segmentFile.bed>           |Your output segmentFile.bed (BED9 format)

**Optional Arguments**

|Flag|Type|Desription|
|----|----|----------|
|-n  --threads      | number                                |number of processors, default = 1

An example of the command is therefore:

    $ FStitch segment -b </path/to/bedGraphFile> -s (+/-) -p </path/to/Parameters.hmminfo> -o </path/to/segmentFile.bed>

This will produce a file called segmentFile.bed, and can be imported into any genome browser)

Note you can use your parameter out file from FStitch train (i.e. Parameters.hmminfo) to segment other datasets in the same series of samples from an experiment. In fact, using the same parameter out file will gurantee consistency and comparibility across datasets, so this is encouraged. 

*That said*, remember that FStitch attempts to discern "noise" from "signal". Therefore, if you have a sample that has low complexity (i.e. it is difficult to distinguish signal from noise), you may either want to train using other samples, perform deeper sequencing on this sample (i.e. obtain a technical replicate -- this may be especially necessary depending on your downstream analysis interests), or train this sample separately (reduces strength of a statistical comparison between samples). As such, sample of comperable complexity and read depth will typically behave the best in FStitch.

## Cite
If you find the Fast Read Stitcher program useful for your research please cite:

Joseph Azofeifa, Mary A. Allen, Manuel Lladser, and Robin Dowell. 2014. __FStitch: a fast and simple algorithm for detecting nascent RNA transcripts.__ In Proceedings of the 5th ACM Conference on Bioinformatics, Computational Biology, and Health Informatics (BCB '14). ACM, New York, NY, USA, 174-183. 

##References 
1. Joseph Azofeifa, Mary A. Allen, Manuel Lladser, and Robin Dowell. 2014. __FStitch: a fast and simple algorithm for detecting nascent RNA transcripts.__ In Proceedings of the 5th ACM Conference on Bioinformatics, Computational Biology, and Health Informatics (BCB '14). ACM, New York, NY, USA, 174-183. DOI=10.1145/2649387.2649427 http://doi.acm.org/10.1145/2649387.2649427   
2. http://genome.ucsc.edu/goldenpath/help/bedgraph.html
3. http://bedtools.readthedocs.org/en/latest/
4. http://openmp.org/wp/

