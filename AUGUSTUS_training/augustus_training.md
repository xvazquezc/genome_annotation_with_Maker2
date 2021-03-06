AUGUSTUS training protocol
================
Xabier Vázquez-Campos
2018-03-06

-   [Pre-requisites](#pre-requisites)
    -   [Software](#software)
        -   [Important](#important)
    -   [Files](#files)
-   [Directory structure](#directory-structure)
-   [Variables](#variables)
-   [Modifications to augustus training script files](#modifications-to-augustus-training-script-files)
    -   [autoAug.pl](#autoaug.pl)
    -   [autoAugTrain.pl](#autoaugtrain.pl)
        -   [Changes](#changes)
    -   [autoAugPred.pl](#autoaugpred.pl)
        -   [Changes](#changes-1)

Pre-requisites
--------------

### Software

-   AUGUSTUS, `conda install augustus`
-   BLAT, or even better, [PBLAT](http://icebert.github.io/pblat/)
-   pslcDNAfilter, `conda install ucsc-pslcdnafilter`

#### Important

`autoAugTrain.pl` does not parse `augustus` properly. To solve this just edit any instance calling the binary.

### Files

-   Genome assembly
-   RNA-seq, assembled (alternatively protein data)
-   Training file, output from `zff2genbank`

Directory structure
-------------------

``` bash
# \
# |-- ${MYGENOME}
# |   |-- maker
# |       |-- augtrain1
# |           |-- cdna
# |           |-- hints
# |           |-- seq
# |           |-- training
```

Variables
---------

``` bash
AUG_TRAIN=~/Desktop/SBI_projects/augustus_training
MYSP=pugra_adv_mpi_v1
GENOME=/home/xabi/Desktop/SBI_projects/Chongmei/pugra_assembly_rename.fasta
EST=/home/xabi/Desktop/SBI_projects/Chongmei/pugra_rna_rename.fasta
TRAIN=~/Desktop/SBI_projects/augustus_training/pugra.train1.gb
```

``` bash
cd ${AUG_TRAIN}

mkdir -p seq
mkdir -p cdna
mkdir -p hints

cd ${AUG_TRAIN}/seq

ln -s ${GENOME} genome.fa
cleanDOSfasta.pl genome.fa > genome_clean.fasta

summarizeACGTcontent.pl genome_clean.fasta > genome.summary

cat genome.summary | grep "bases." | perl -pe 's/(\d+)\sbases.\s+(\S*) BASE.*/$2\tassembly\tcontig\t1\t$1\t.\t.\t.\tContig $2/' > contigs.gff

cd ${AUG_TRAIN}/cdna
pblat -noHead -threads=${PBS_NUM_PPN} -minIdentity=80 -maxIntron=3000 ${AUG_TRAIN}/seq/genome_clean.fasta ${EST} cdna.psl 1>blat.stdout 2>blat.stderr
# original -maxIntron value was 100000, which is absurdly high for fungal genomes. Set to 3000 as for HISAT2

pslCDnaFilter -minId=0.9 -localNearBest=0.005 -ignoreNs -bestOverlap cdna.psl cdna.f.psl 1>pslCDnaFilter.stdout 2>pslCDnaFilter.stderr

cd ${AUG_TRAIN}/hints

blat2hints.pl --in=${AUG_TRAIN}/cdna/cdna.f.psl --out=hints.gff --minintronlen=35 --trunkSS 1>blat2hints.stdout 2>blat2hints.stderr


####### Step 1: Training AUGUSTUS (no UTR models) #######
echo -e "####### Step 1: Training AUGUSTUS (no UTR models) #######"
cd ${AUG_TRAIN}
autoAugTrain.pl --trainingset=${TRAIN} --species=${MYSP} --genome=${AUG_TRAIN}/seq/genome_clean.fasta --workingdir=${AUG_TRAIN} -v -v -v --optrounds=1 --cpus=${PBS_NUM_PPN}


####### Step 2: Preparing scripts for AUGUSTUS without hints and UTR #######
echo -e "####### Step 2: Preparing scripts for AUGUSTUS without hints and UTR #######"
autoAugPred.pl --genome=${AUG_TRAIN}/seq/genome_clean.fasta --species=${MYSP} --workingdir=${AUG_TRAIN} -v -v -v --singleCPU


####### Step 3: Continue to predict genome structure with AUGUSTUS without hints, no UTR #######
autoAugPred.pl --species=${MYSP} --genome=${AUG_TRAIN}/seq/genome_clean.fasta --continue --workingdir=${AUG_TRAIN} -v -v -v --singleCPU


####### Step 4: Preparing scripts for AUGUSTUS with hints, without UTR #######
autoAugPred.pl --genome=${AUG_TRAIN}/seq/genome_clean.fasta --species=${MYSP} --workingdir=${AUG_TRAIN} -v -v -v --hints=${AUG_TRAIN}/hints/hints.gff --singleCPU
# if done in the cluster, do cd to the ${AUG_TRAIN}/autoAugPred_hints/shells folder, before running. autoAugPred.pl passes the hints file with a relative path to each of the parallel files

####### Step 5: Continue to predict genome structure with AUGUSTUS with hints, no UTR #######
autoAugPred.pl --species=${MYSP} --genome=${AUG_TRAIN}/seq/genome_clean.fasta --continue --workingdir=${AUG_TRAIN} -v -v -v --hints=${AUG_TRAIN}/hints/hints.gff --singleCPU


####### Step 6: Training AUGUSTUS with UTR #######
autoAugTrain.pl --cpus=4 --genome=${AUG_TRAIN}/seq/genome_clean.fasta --species=${MYSP} --utr --estali=${AUG_TRAIN}/cdna/cdna.f.psl --aug=${AUG_TRAIN}/autoAugPred_hints/predictions/augustus.gff --workingdir=${AUG_TRAIN} -v -v -v --opt=1 --useexisting --cpus=${PBS_NUM_PPN}


####### Step 7: Preparing scripts for AUGUSTUS without hints, with UTR #######
autoAugPred.pl --species=${MYSP} --genome=${AUG_TRAIN}/seq/genome_clean.fasta --continue --utr --workingdir=${AUG_TRAIN} -v -v -v --singleCPU

####### Step 8: Continue to predict genome structure with AUGUSTUS without hints, with UTR #######
autoAugPred.pl --species=${MYSP} --genome=${AUG_TRAIN}/seq/genome_clean.fasta --continue --utr --workingdir=${AUG_TRAIN} -v -v -v --singleCPU

####### Step 9: Preparing scripts for AUGUSTUS with hints and UTR #######
autoAugPred.pl --genome=${AUG_TRAIN}/seq/genome_clean.fasta --species=${MYSP} --utr --workingdir=${AUG_TRAIN} -v -v -v --hints=${AUG_TRAIN}/hints/hints.gff --singleCPU
# if done in the cluster, do cd to the ${AUG_TRAIN}/autoAugPred_hints/shells folder, before running. autoAugPred.pl passes the hints file with a relative path to each of the parallel files

####### Step 10: Continue to predict genome structure with AUGUSTUS with hints and UTR #######
autoAugPred.pl --species=${MYSP} --genome=${AUG_TRAIN}/seq/genome_clean.fasta --continue --utr --workingdir=${AUG_TRAIN} -v -v -v --hints=${AUG_TRAIN}/hints/hints.gff --singleCPU
```

Modifications to augustus training script files
-----------------------------------------------

### autoAug.pl

Uses BLAT to create the hints file. BLAT runs in single cpu and it's pretty slow. We substitute BLAT for PBLAT and we add the `--cpus` option. This will be passed into PBLAT, obviously, but also to `autoAugTrain.pl`, which with additional mods, will pass the option into `optimize_augustus.pl`.
\#\#\#\# Changes Add at line 55:

``` bash
my $cpus=1;                           # number of cpus for pblat and optimize_augustus.pl
```

Add line 108:

``` bash
--cpus=number                       Number of CPUs to use in PBLAT and optimize_augustus.pl.
```

Line 135, within `GetOptions`:

``` bash
'cpus=i' => \$cpus
```

Change line 674:

``` bash
$cmdString="pblat -threads=$cpus -noHead  -minIdentity=80 -maxIntron=$maxIntronLen ../seq/genome_clean.fa cdna.fa cdna.psl 1>blat.stdout 2>blat.stderr"; 
```

At line 756, add the `--cpus` option:

``` bash
$perlCmdString="perl $scriptPath/autoAugTrain.pl --cpus=$cpus -t=$trainingset -s=$species $useexistingopt -g=$genome_clean -w=$rootDir $verboseString --opt=$optrounds";
```

Change lines 917-921 to this (basically add support for `--cpus`)

``` bash
    if(-d $rootDir){
      $perlCmdString="perl $scriptPath/autoAugTrain.pl --cpus=$cpus -g=$genome_clean -s=$species --utr -e=$estali $augString -w=$rootDir $verboseString --opt=$optrounds --useexisting";
    }else{
      $perlCmdString="perl $scriptPath/autoAugTrain.pl --cpus=$cpus -g=$genome_clean -s=$species --utr -e=$estali $augString -w=$rootDir $verboseString --opt=$optrounds $useexistingopt";
    }
```

Change `blat` for `pblat` in the installation checks (line 1084):

``` bash
    if (system("which pblat > /dev/null") != 0){
```

### autoAugTrain.pl

-   Bug in the generation of UTR evidence
-   Doesn't provide a cpu option to be passed into `optimizeAugustus.pl`

#### Changes

Add at L45:

``` bash
my $cpus=1;                # number of cpus for optimize_augustus.pl
```

Add L85:

``` bash
$usage.="--cpus=number                  Number of CPUs to use in optimize_augustus.pl.\n";
```

Add at L93 (within GetOptions):

``` bash
       'cpus=i' => \$cpus,
```

Modify L378-382 to pass through the `--cpus` option to `optimizeAugustus.pl`:

``` bash
    if($t_b_o==0){
        $cmdString="perl $string --cpus=$cpus --rounds=$optrounds --species=$species $workDir/training/training.gb.train.test --metapars=$configDir/$metaName > optimize.out";
    } else{
        $cmdString="perl $string --cpus=$cpus --rounds=$optrounds --species=$species $workDir/training/training.gb.train.test --onlytrain=$workDir/training/training.gb.onlytrain --metapars=$configDir/$metaName > optimize.out";
    }
```

At L559-569, substitute the lines between the `open` and `close` statements. This is what causes the UTR training to fail!!:

``` bash
    open(TR, "tr.lst") or die ("Can not open tr.lst!\n");
    open(BOTH, "> bothutr.lst");
    my $Fld1;
    my $prev;
    while(<TR>){
        ($Fld1) = split('\t', $_, -1);
        if ($Fld1 eq $prev) {   #???
        print BOTH "$prev\n";
        }
        $prev = $Fld1;
    }
    close(TR);
    close(BOTH);
```

At L716, add the last `--cpus` option to `optimizeAugustus.pl`:

``` bash
    $perlCmdString="perl $string --cpus=$cpus --rounds=$optrounds --species=$species --trainOnlyUtr=1 --onlytrain=onlytrain.gb  --metapars=$configDir/$metaUtrName train.gb --UTR=on > optimize.utr.out";
```

### autoAugPred.pl

Has SGE paths hard-coded, need to be changed to be able to use Katana (Torque/PBS) capabilities and not having to rely on the `--singleCPU` mode.

#### Changes

Change L46-47:

``` bash
my $clusterENVDefs = "export SGE_ROOT=/opt/torque";  # commands to define environment variables on cluster
my $SGEqstatPath = "/opt/torque/bin/";    # path to executables on Sun Grid Engine (qsub, qstat)
```

L272 has given problems before, change to this:

``` bash
        $opt_string .= " --hintsfile=$workDir/../hints/hints.gff --extrinsicCfgFile=$extrinsiccfg"; 
```

If you use AUGUSTUS installed from conda you need to change L275:

``` bash
        my $aug="$AUGUSTUS_CONFIG_PATH/../bin/augustus";
```

For Katana, we need to remove the `-cwd` option from the `qsub` command at L296:

``` bash
        print SH '    qsub "aug$i"'."\n";  #just for sun grid engine
```

In some cases, the default memory for Katana jobs caused crashes, so I changed the L483 to this. Probably a bit overkill to request 4gb, but that's a safe choice:

``` bash
    my $string_2="\"cd $workDir; $clusterENVDefs;${SGEqstatPath}qsub -l vmem=4gb -e aug$j.err ";
```

The job id in Katana follows this structure: `6481541.katana.science.unsw.edu.au`, but in the output of `qstat`, we usually only see `6481541.katana`, so we need to add some lines so the script is able to find it when it runs the perl subroutine `waitJob`, which uses `grep` in the `qstat` output to evaluate if the jobs have already finished or not. This applies to L488-497.
In addition, `autoAugPred.pl` expects `/^Your job (\d+)/` in the `stdout` of `qsub` and uses it to grab the job ID, but Katana throws the job ID directy. This is at L501.
This chunk covers L486-508:

``` bash
    print "3 $cmd\n" if ($verbose>=3);
    system("$cmd")==0 or die ("failed to execute: $!\n");
    open(J, "qsub$j.stdout") or die ("Error: could not open a qsub stdout file to fix your job number!\n");
    open(JFIX, "> qsub$j.stdout_fix");
    while(my $jobid = <J>){
        chomp $jobid;
        my @col = split(/\./, $jobid), "\n";
        print JFIX join(".", $col[0], $col[1]), "\n";
    }
    close(J);
    close(JFIX);
    $cmd="cat qsub$j.stdout_fix";
    my $output = `$cmd 2>/dev/null`;
    print "3 $output" if ($verbose>=3);
        # exception (???)
    if($output=~/^(\d+)/){
        print TEMP $1."\n";
    }
    $j++;
    sleep(1);
  }
  close TEMP;
  waitJob($workDir, "your_job.tmp");
```
