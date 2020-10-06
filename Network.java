/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


package sim;

import cern.jet.random.Poisson;
import cern.jet.random.Normal;
import cern.jet.random.engine.DRand;
import cern.jet.random.engine.RandomEngine;
import com.jmatio.io.MatFileReader;
import com.jmatio.types.MLDouble;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.Math.exp;
import static java.lang.Math.random;
import static java.lang.Math.sqrt;
import static java.lang.Math.log;
import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.pow;
import static java.lang.Math.ceil;
import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.Math.signum;
import static java.lang.Math.abs;
/**
 *
 * @author Logan
 */
public class Network {
    /*===Simulation time parameters===*/
    public double t = 0; //ms
    public double dt = 0.1; //ms

    /*===Basic network layout parameters===*/
    public int numHCRows,numHC,numExcitatoryPerHC,numInhibitoryPerHC; //HC = "hypercolumn"; model consists of a square array of hypercolumns
    public int numExcitatory,numInhibitory;
    public int numLGNBlockRows,numLGNBlocks; //LGN is organized into a square array of blocks
    public double numLGNPerRowPerBlock,numLGNPerColPerBlock;
    public int numLGNPerBlock,numLGN;
    public int HCWidth,patchWidth;
    public MyPoint[] locExc,locInh,locLGNON,locLGNOFF;

    /*===Connectivity parameters and matrices===*/
    public int[] numEPostExc,numIPostExc,numEPostInh,numIPostInh,numEPreExc,numIPreExc,numEPreInh,numIPreInh; //numEPostExc[i] = number of E cells postsynaptic to the ith E cell; numEPreExc[i] = number of E cells presynaptic to ith E cell
    public int[][] cEExc,cIExc,cEInh,cIInh,cPreEExc,cPreIExc,cPreEInh,cPreIInh; //cIExc[i][j] = index of jth postsynaptic I cell to the E cell at index i; cPreIExc[i][j] = index of jth presynaptic I cell to the E cell at index i

    /*===Connection strengths===*/
    public double sLGN,sLGNInh,sEE,sIE,sEI,sII, justSEE;
    public double sModSTD = 0.1;
    public double[] sModSeedExc, sModSeedInh, sModSeedLGNON, sModSeedLGNOFF;

    /*===Layer 6 recording parameters===*/
    public double[] fbSpikeTimes;
    public int[] fbSpikes;
    public int numFBSpikes, numFB;
    public int[][] excPostFB, inhPostFB;
    public int[] numExcPostFB, numInhPostFB;
    public double[][] locFB;
    public double fbClock = 400;
    public double fbStart = 400;
    public double fbEnd = 900;
    public String fbInputFile;
    public int currentFBSpike = 0;
    public double[] fbExtraProb;
    public double[] fbInitProb;
    public double fbISD = 85.0;
    public double fbESD = 110.0;
    public double fbIPeak = .6;
    public double pHigh = 2.0;
    public double pLow = 0.0;
    public double[] fbPValue;
    public double[] pValues = { 2.0, 1.6, 1.6, 1.2, 1.2, 0.8, 0.8, 1.2, 1.2, 1.6, 1.6, 2.0, 2.0 };//{ 2.0, 1.5, 1.5, 1.0, 1.0, 0.5, 0.5, 1.0, 1.0, 1.5, 1.5, 2.0, 2.0 };

    /*===Default Network parameters===*/
    public double gLExc = 50.0/1000.0;
    public double gLInh = gLExc / 0.75;
    public double eRev = 4.67;
    public double iRev = -0.67;
    public double refractoryE = 2.0; //ms
    public double refractoryI = 1.0;
    public double tauAMPA = 1;
    public double tauGABA = 1.67;
    public double tauNMDA1 = 80.0;
    public double tauNMDA2 = 2.0;
    public double baseLGNRate = 20.0 / 1000.0; //spikes-per-ms
    public double eLenScale = 200.0 / sqrt(2);
    public double iLenScale = 125.0 / sqrt(2);

    /*===Neuron state variables===*/
    public double[] vExc,vInh,sleepExc,sleepInh; //E and I membrane potentials, and refractory period clocks
    public double[][] gLGNExc,gLGNInh,gAMPAExc,gAMPAInh,gGABAExc,gGABAInh,gNMDAExc,gNMDAInh; //conductances

    /*===LGN -- Cortex connectivity===*/
    public int[][] cELGNON; //cELGNON[i][j] = index of jth E cell postsynaptic to ith LGN ON cell
    public int[][] cILGNON;
    public int[][] cELGNOFF;
    public int[][] cILGNOFF;
    public int[] nELGNON; //nELGNON[i] = number of postsynaptic E connections from ith LGN ON cell
    public int[] nILGNON;
    public int[] nELGNOFF;
    public int[] nILGNOFF;
    public int[] numInputsExc; //numInputsExc[i] = number of LGN ON and OFF inputs to excitatory cell i
    public int[] numInputsInh;
    public String LGNInputFile; //LGN connectivity data is placed in a MATLAB .mat file, the name of which is stored here

    /*===Temporary addin for low density experiments (currently supercedes connectivityQQ' variables)===*/
    public double[] plainsEE; //"plains refers to the neurons away from the boundary of the simulation"
    public double[] plainsEI; 
    public double[] plainsIE;
    public double[] plainsII;
    public double[] bndryEE;
    public double[] bndryEI;
    public double[] bndryIE;
    public double[] bndryII;

    /*===egalitarian feedback multipliers===*/
    public double[] fbEgalMod;

    /*== LGN parameters ==*/
    public double leakLGN = 1.0/10.0;
    public double[] phaseON, phaseOFF, vLGNON,vLGNOFF;
    public double spFreqMod;

    /*===LGN Drift grating parameters===*/
    public double I_0 = 100.0 / 1000.0;
    public double shotNoise = 0.08;
    public double noiseRate = 1.0;
    public double threshLGN = 1.15;
    public double eps = 2.0;
    public double w = 4.0/1000.0;
    public double kFreq = .025;
    public MyPoint kOrient = new MyPoint(1,0); //MyPoint is just an extension of the point2D class
    public double kAngle = 0.0;

    /*===Extra Drive parameters (NOTE: extra drive conductance will be a part of LGN conductance)===*/
    public double rateStrongAmbE;
    public double rateStrongAmbI;
    public double rateWeakAmbE;
    public double rateWeakAmbI;
    public double sStrongAmbE = sEE;
    public double sStrongAmbI = sIE;
    public double sWeakAmb = 0.01;
//    public double rateTotalE;
//    public double rateTotalI;
//    public double strongProbE;
//    public double strongProbI;

    /*===Feedback parameters===*/
    public double fbPref = 0;
    public double fbOrth = 0;
    public double fbAvg;
    public double fbMod;


    /*===Quick stats collecting mechanism===*/
    public int statSkipInSteps = 20;
    public int statStep = 0; //will be incremented mod statStkipsInSteps every step
    public double patchOneFRExc = 0;
    public double patchOneFRExcTemp = 0;
    public double patchOneFRInh = 0;
    public double patchOneFRInhTemp = 0;
    public double patchTwoFRExc = 0;
    public double patchTwoFRExcTemp = 0;
    public double patchTwoFRInh = 0;
    public double patchTwoFRInhTemp = 0;
    public int[] membershipExc;
    public int[] membershipInh;

    /*===Random number generators===*/
    public RandomEngine engine = new DRand(new java.util.Date());
    public Poisson poisson = new Poisson(1.0, engine);
    public Normal normal = new Normal(0.0,1.0, engine);

    /*===Stats===*/
    public int eSpikes = 0;
    public int iSpikes = 0;
    public int numSteps = 0;
    public int[] spikesNowExc;
    public int[] spikesNowInh;
    public double[] spikeTimesNowExc;
    public double[] spikeTimesNowInh;
    public int numSpikesNowExc;
    public int numSpikesNowInh;

    /*===I randomized threshold===*/
    public double[] threshI;
    public double[] threshE;

    /*===I input modifier===*/
    public double modIILow = .9;
    public double modIIHigh = 1.1;
    public double modIIRange;

    /*===Fast spike initial conductance for strength 1===*/
    public double[] fastSpikeInitial = {0.041193481496372,0.037000133080574,0.033233652467581,0.029850586048726};//1.5sec speedup////{0.055017482305651,0.032944600183025,0.019727305498817,0.011812757783723};//1 sec speedup//{0.047552705667974,0.035593342565849,0.026641723477432,0.019941409788497};//1.25 sec speedup////nospeedup:{ 0.093999733184618, 0.005628726538001, 0.000337049493293, 0.000020182604389};//{0.041177,0.036986,0.033221,0.029839}; //exp,lin-exp,quad-exp,cub-exp

    /*===E->E spike delay cylinder (a data structure for quickly implementing E hits with randomized delays)===*/
    public double spikeDelayMs = 3.0;
    public int spikeDelaySteps = (int) ceil(spikeDelayMs / dt);
    public double[][] delayCylinder;
    public int spikeDelayIndex = 0;
    public int lgnSpikeDelaySteps;
    public double sLGNOverSEEFactor;

    /*=== for degugging ===*/
    public int debugCount=0;

    /*===Constructor===*/
    public Network(int numHCRows, int numExcitatoryPerHC, int numInhibitoryPerHC, int numLGNBlockRows, double numLGNPerRowPerBlock, double numLGNPerColPerBlock){
        //Network is not built in constructor--this is done with buildNetwork() and allocateConnectivity().  Many parameters are set from MATLAB after this is called, BEFORE connectivity is implemented.
        this.numHCRows = numHCRows;
        numHC = numHCRows * numHCRows;
        this.numExcitatoryPerHC = numExcitatoryPerHC;
        this.numInhibitoryPerHC = numInhibitoryPerHC;
        numExcitatory = numExcitatoryPerHC*numHC;
        numInhibitory = numInhibitoryPerHC*numHC;

        this.numLGNBlockRows = numLGNBlockRows;
        numLGNBlocks = numLGNBlockRows * numLGNBlockRows;
        this.numLGNPerRowPerBlock = numLGNPerRowPerBlock;
        this.numLGNPerColPerBlock = numLGNPerColPerBlock;
        this.numLGNPerBlock = (int) Math.round(numLGNPerRowPerBlock * numLGNPerColPerBlock);
        numLGN = numLGNBlocks * numLGNPerBlock;

        HCWidth = 500; //in micrometers
        patchWidth = HCWidth * numHCRows;

        /*===Network state===*/
        vExc = new double[numExcitatory];
        vInh = new double[numInhibitory];
        sleepExc = new double[numExcitatory];
        sleepInh = new double[numInhibitory];

        /*==LGN state==*/
        vLGNON = new double[numLGN];
        vLGNOFF= new double[numLGN];
        phaseON= new double[numLGN];
        phaseOFF=new double[numLGN];


        /*===Network state helpers===*/
        gLGNExc = new double[numExcitatory][4];
        gLGNInh = new double[numInhibitory][4];
        gAMPAExc = new double[numExcitatory][4];
        gAMPAInh = new double[numInhibitory][4];
        gGABAExc = new double[numExcitatory][4];
        gGABAInh = new double[numInhibitory][4];
        gNMDAExc = new double[numExcitatory][3];
        gNMDAInh = new double[numInhibitory][3];

        /*=== Place neuron locations ===*/
        int ePopWidth = (int) ceil(sqrt(numExcitatory)); //in neurons
        int iPopWidth = (int) ceil(sqrt(numInhibitory)); //in neurons
        locExc = new MyPoint[numExcitatory];
        locInh = new MyPoint[numInhibitory];
        for (int i=0;i<numExcitatory;i++){
            locExc[i] = new MyPoint((double)(i % ePopWidth) * (patchWidth/(double)ePopWidth),
                                    (double)(i / ePopWidth) * (patchWidth/(double)ePopWidth));
        }
        for(int i = 0; i < numInhibitory; i++){
            locInh[i] = new MyPoint((double)(i % iPopWidth) * (patchWidth/(double)iPopWidth),
                                    (double)(i / iPopWidth) * (patchWidth/(double)iPopWidth));
        }
        /*=== Neurons placed ===*/

        membershipExc = new int[numExcitatory];
        membershipInh = new int[numInhibitory];

        spikesNowExc = new int[numExcitatory];
        spikesNowInh = new int[numInhibitory];
        spikeTimesNowExc = new double[numExcitatory];
        spikeTimesNowInh = new double[numInhibitory];

        threshI = new double[numInhibitory];
        threshE = new double[numExcitatory];

        sModSeedExc = new double[numExcitatory];
        sModSeedInh = new double[numInhibitory];
        sModSeedLGNON = new double[numLGN];
        sModSeedLGNOFF= new double[numLGN];

        delayCylinder = new double[numExcitatory][spikeDelaySteps];

        numInputsExc = new int[numExcitatory];
        numInputsInh = new int[numInhibitory];        

        nELGNON = new int[numLGN]; //read: "number connections to E from LGN ON"
        nILGNON = new int[numLGN]; //read: "number connections to E from LGN ON"
        nELGNOFF = new int[numLGN]; //read: "number connections to E from LGN ON"
        nILGNOFF = new int[numLGN]; //read: "number connections to E from LGN ON"
    }

	/* allocate space to place the connectivity matrix (and do some other initializations) */
    public void allocateConnectivity(double connectivityEE,double connectivityIE,double connectivityEI, double connectivityII, double[] egalEHelper, double[] egalIHelper, String LGNInputFile, String fbInputFile){
        this.LGNInputFile = LGNInputFile;
        this.fbInputFile = fbInputFile;

        /*===Calculations to improve efficient memory allocation for connectivity matrix===*/
        double egalEBump = egalEHelper[0];
        double egalIBump = egalIHelper[0];
        double EDensity = (double)numExcitatoryPerHC / (double)(HCWidth*HCWidth);
        double IDensity = (double)numInhibitoryPerHC / (double)(HCWidth*HCWidth);
        double EDiscArea = PI * pow(2 * eLenScale,2);
        double IDiscArea = PI * pow(2 * iLenScale,2);
        double numEInEDisc = EDiscArea * EDensity;
        double numIInEDisc = EDiscArea * IDensity;
        double numEInIDisc = IDiscArea * EDensity;
        double numIInIDisc = IDiscArea * IDensity;
        int maxEToEOutDegree = (int)floor((connectivityEE+egalEBump) * numEInEDisc);
        int maxEToIOutDegree = (int)floor((connectivityIE+egalEBump) * numIInEDisc);
        int maxIToEOutDegree = (int)floor((connectivityEI+egalIBump) * numEInIDisc);
        int maxIToIOutDegree = (int)floor((connectivityII+egalIBump) * numIInIDisc);
        int maxEToEInDegree = maxEToEOutDegree;
        int maxEToIInDegree = (int)floor((connectivityIE+egalEBump) * numEInEDisc);
        int maxIToEInDegree = (int)floor((connectivityEI+egalIBump) * numIInIDisc);
        int maxIToIInDegree = maxIToIOutDegree;

        System.out.println(maxEToEOutDegree);
        System.out.println(maxEToIOutDegree);
        System.out.println(maxIToEOutDegree);
        System.out.println(maxIToIOutDegree);
        System.out.println(maxEToEInDegree);
        System.out.println(maxEToIInDegree);
        System.out.println(maxIToEInDegree);
        System.out.println(maxIToIInDegree);

        /*===Network Connectivity===*/
        numEPostExc = new int[numExcitatory]; //num E cells postsynaptic to excitatory [array element]
        numIPostExc = new int[numExcitatory];
        numEPostInh = new int[numInhibitory];
        numIPostInh = new int[numInhibitory];
        numEPreExc = new int[numExcitatory];
        numIPreExc = new int[numExcitatory];
        numEPreInh = new int[numInhibitory];
        numIPreInh = new int[numInhibitory];
        cEExc = new int[numExcitatory][maxEToEOutDegree];
        cIExc = new int[numExcitatory][maxEToIOutDegree];
        cEInh = new int[numInhibitory][maxIToEOutDegree];
        cIInh = new int[numInhibitory][maxIToIOutDegree];
        cPreEExc = new int[numExcitatory][maxEToEInDegree];
        cPreIExc = new int[numExcitatory][maxIToEInDegree];
        cPreEInh = new int[numInhibitory][maxEToIInDegree];
        cPreIInh = new int[numInhibitory][maxIToIInDegree];
        locLGNON = new MyPoint[numLGN];
        locLGNOFF = new MyPoint[numLGN];
    }

    public void buildNetwork(){
        /* Initialize all fields */

        /*== Get LGN inputs ==*/
        getLGNInputs();
        /*== Got LGN inputs ==*/

        /*== setup cortical network ==*/
        setupCortexNetwork();
        /*== cortical network setup ==*/

        /*== Initialize I,E thresholds ==*/
        for (int i = 0;i<numInhibitory;i++){
            threshI[i] = 1.0 - 0.2 * random();
        }
        for (int i = 0;i<numExcitatory;i++){
            threshE[i] = 1.0;
        }

        /*== set strength modifiers ==*/
        setSModSeeds();
    }

    public void setupCortexNetwork(){
        /*=== make sure degree numbers are set to 0 ===*/
        for (int i = 0; i < numExcitatory; i++) {
            numEPreExc[i] = 0;
            numEPostExc[i] = 0;
            numIPreExc[i] = 0;
            numIPostExc[i] = 0;
        }
        for (int i = 0; i < numInhibitory; i++) {
            numEPreInh[i] = 0;
            numEPostInh[i] = 0;
            numIPreInh[i] = 0;
            numIPostInh[i] = 0;
        }
        /*=== degrees reset ===*/

        /* calculate expected number presynapic for each pair type*/
        double EDensity = (double)numExcitatoryPerHC / ((double)HCWidth * (double)HCWidth); // num E per square micrometer
        double IDensity = (double)numInhibitoryPerHC / ((double)HCWidth * (double)HCWidth); // num E per square micrometer

        double EDensityPerESD = EDensity * eLenScale * eLenScale;
        double EDensityPerISD = EDensity * iLenScale * iLenScale;
        double IDensityPerESD = IDensity * eLenScale * eLenScale;
        double IDensityPerISD = IDensity * iLenScale * iLenScale;

        double desiredNumPrePeak1EToE = 2.0*PI*(1.0-exp(-2.0)) * EDensityPerESD;
        double desiredNumPrePeak1IToE = 2.0*PI*(1.0-exp(-2.0)) * IDensityPerISD;
        double desiredNumPrePeak1EToI = 2.0*PI*(1.0-exp(-2.0)) * EDensityPerESD;
        double desiredNumPrePeak1IToI = 2.0*PI*(1.0-exp(-2.0)) * IDensityPerISD;

        System.out.println("Typical Pre:");
        System.out.println("E->E: " + (desiredNumPrePeak1EToE*.15));
        System.out.println("I->E: " + (desiredNumPrePeak1IToE*.6));
        System.out.println("E->I: " + (desiredNumPrePeak1EToI*.6));
        System.out.println("I->I: " + (desiredNumPrePeak1IToI*.6));

        /*== E->E and I->E ==*/
        for (int i=0;i<numExcitatory;i++){
            //choosePresynaptic(postIndex,postLocation,sizeOfPreGroup,preGroupLocs, lenScale, expNumPre, cOut, numCOut, cIn, numCIn)
            /*== choose E presynaptic ==*/
            /*== first compute the required bump to connectivity to keep expected num presynaptic the same ==*/
            if (i%3000 == 0) System.out.println("Connecting E neuron " + i);
            choosePresynaptic(i,
                    locExc[i],
                    numExcitatory,
                    locExc,
                    eLenScale,
                    numInputsExc,
                    cEExc,
                    numEPostExc,
                    cPreEExc,
                    numEPreExc,
                    plainsEE,
                    bndryEE,
                    desiredNumPrePeak1EToE);

            /*== choose I presynaptic ==*/
            choosePresynaptic(i,
                    locExc[i],
                    numInhibitory,
                    locInh,
                    iLenScale,
                    numInputsExc,
                    cEInh,
                    numEPostInh,
                    cPreIExc,
                    numIPreExc,
                    plainsEI,
                    bndryEI,
                    desiredNumPrePeak1IToE);
        }

        /*== redraw if outside 10% and 90% range ==*/
        System.out.println("redrawing E->E and I->E outliers...");
        for (int numLGNPre=0; numLGNPre<9; numLGNPre++){
            System.out.println("numLGN=" + numLGNPre);
            //find all E cells with numLGNPre LGN inputs
            int[] picks = new int[numExcitatory];
            int numPicks = 0;
            for (int i=0;i<numExcitatory;i++){
                if (numInputsExc[i]==numLGNPre){
                    picks[numPicks]=i;
                    numPicks++;
                }
            }
            //found all E cells with numLGNPre LGN inputs

            //get list of all the number presynaptic
            Pair[] preSynEPicks = new Pair[numPicks];
            Pair[] preSynIPicks = new Pair[numPicks];
            for (int i=0;i<numPicks;i++){
                preSynEPicks[i] = new Pair(numEPreExc[picks[i]], picks[i]);
                preSynIPicks[i] = new Pair(numIPreExc[picks[i]], picks[i]);
            }

            //sort the list of number presynaptic
            Arrays.sort(preSynEPicks);
            Arrays.sort(preSynIPicks);

            //redraw all <10%
            for (int i=0;i<(int)floor(0.1*(double)numPicks);i++){
                deletePresynaptic(preSynEPicks[i].index,
                        cEExc,
                        numEPostExc,
                        cPreEExc,
                        numEPreExc);
                choosePresynaptic(preSynEPicks[i].index,
                    locExc[preSynEPicks[i].index],
                    numExcitatory,
                    locExc,
                    eLenScale,
                    numInputsExc,
                    cEExc,
                    numEPostExc,
                    cPreEExc,
                    numEPreExc,
                    plainsEE,
                    bndryEE,
                    desiredNumPrePeak1EToE);

                deletePresynaptic(preSynIPicks[i].index,
                        cEInh,
                        numEPostInh,
                        cPreIExc,
                        numIPreExc);
                choosePresynaptic(preSynIPicks[i].index,
                    locExc[preSynIPicks[i].index],
                    numInhibitory,
                    locInh,
                    iLenScale,
                    numInputsExc,
                    cEInh,
                    numEPostInh,
                    cPreIExc,
                    numIPreExc,
                    plainsEI,
                    bndryEI,
                    desiredNumPrePeak1IToE);
            }

            //redraw all >90%
            for (int i=(int)floor(0.9*(double)numPicks);i<numPicks;i++){
                deletePresynaptic(preSynEPicks[i].index,
                        cEExc,
                        numEPostExc,
                        cPreEExc,
                        numEPreExc);
                choosePresynaptic(preSynEPicks[i].index,
                    locExc[preSynEPicks[i].index],
                    numExcitatory,
                    locExc,
                    eLenScale,
                    numInputsExc,
                    cEExc,
                    numEPostExc,
                    cPreEExc,
                    numEPreExc,
                    plainsEE,
                    bndryEE,
                    desiredNumPrePeak1EToE);

                deletePresynaptic(preSynIPicks[i].index,
                        cEInh,
                        numEPostInh,
                        cPreIExc,
                        numIPreExc);
                choosePresynaptic(preSynIPicks[i].index,
                    locExc[preSynIPicks[i].index],
                    numInhibitory,
                    locInh,
                    iLenScale,
                    numInputsExc,
                    cEInh,
                    numEPostInh,
                    cPreIExc,
                    numIPreExc,
                    plainsEI,
                    bndryEI,
                    desiredNumPrePeak1IToE);
            }
        }
        System.out.println("done.");

        /*== E->I and I->I ==*/
        for (int i=0;i<numInhibitory;i++){
            if (i%3000 == 0) System.out.println("Connecting I neuron " + i);
            /*== choose E presynaptic ==*/
            choosePresynaptic(i,
                    locInh[i],
                    numExcitatory,
                    locExc,
                    eLenScale,
                    numInputsInh,
                    cIExc,
                    numIPostExc,
                    cPreEInh,
                    numEPreInh,
                    plainsIE,
                    bndryIE,
                    desiredNumPrePeak1EToI);

            /*== choose I presynaptic ==*/
            choosePresynaptic(i,
                    locInh[i],
                    numInhibitory,
                    locInh,
                    iLenScale,
                    numInputsInh,
                    cIInh,
                    numIPostInh,
                    cPreIInh,
                    numIPreInh,
                    plainsII,
                    bndryII,
                    desiredNumPrePeak1IToI);
        }
        /*== Chose all cortical connections ==*/

        /*== Choose FB connections ==*/
        loadFB();
        /*== FB loaded ==*/
    }

    /*
    public final int[] pickKOfN(int k, int n){
        int numToPick = k;
        int numRemaining = n;
        int[] picks = new int[k];
        int i = 0;
        while (numToPick > 0){
            double p = (double)numToPick / (double)numRemaining;
            if (random() > p){
                picks[(int)numToPick-1] = i;
                numToPick--;
            }
            numRemaining--;
        }
        return picks;
    }
     * */

    public final double quickNormInv(double alpha){
        return (10/log(41)) * log(1 - log(-log(alpha)/log(2))/log(22));
    }

    public final double modifier(double seedi, double seedj){
        return 1.0 + sModSTD * quickNormInv((seedi * 1024.0 + seedj) % 1.0);
    }

    public void deletePresynaptic(int post,
            int[][] cPostFromPre,
            int[] numPostFromPre,
            int[][] cPreToPost,
            int[] numPreToPost){
        for (int i=0; i<numPreToPost[post]; i++){
            int pre = cPreToPost[post][i];

            ////delete connection

            //find post in pre's connectivity list
            int postsPlace=0;
            while(cPostFromPre[pre][postsPlace]!=post){
                postsPlace++;
            }

            //shift all others in list left 1
            for (int j=postsPlace;j<numPostFromPre[pre]-1;j++){
                cPostFromPre[pre][j]=cPostFromPre[pre][j+1];
            }
            numPostFromPre[pre]--;

            ////deleted
        }

        numPreToPost[post]=0;
    }

    public final void choosePresynaptic(int post,
            MyPoint locPost,
            int sizePreGroup,
            MyPoint[] preLocs,
            double lenScale,
            int[] numInputs,
            int[][] cPostFromPre,
            int[] numPostFromPre,
            int[][] cPreToPost,
            int[] numPreToPost,
            double[] egalPlains,
            double[] egalBndry,
            double desiredNumPrePeak1){
        /*== Choose presynaptic connections ==*/
        double s = (double)numInputs[post] / 8.0;
        double distToBndryX = min(locPost.x,patchWidth-locPost.x);
        double distToBndryY = min(locPost.y,patchWidth-locPost.y);
        double distToBndry = min(distToBndryX,distToBndryY);

        double connectivityWithEgal; //to be calculated in following if-else statement

        if (distToBndry < 1.8*lenScale){
            /* then we need to numerically adjust the peak probability to get the expected number presynaptic to be the same */

            /* find all neurons within 2*lenScale */
            double[] allR = new double[3500];
            int numR = 0;
            for (int i=1; i<sizePreGroup;i++){
                double sds = preLocs[i].distance(locPost) / lenScale;
                if (sds < 2.0){
                    allR[numR] = sds;
                    numR++;
                }
            }

            /* map f(x)=e^(-0.5*x^2) over that */
            double[] allP = new double[3500];
            int numP = numR;
            for (int i=1; i<numR; i++){
                allP[i] = exp(-0.5*allR[i]*allR[i]);
            }

            /* use bisection method to find correct connectivity */
            double CLow = egalPlains[numInputs[post]];
            double CHigh = 7.39;
            double expNumPreLow = 0;
            double expNumPreHigh = numR;
            for (int i=1; i<numP; i++){
                expNumPreLow += allP[i];
            }
            expNumPreLow *= CLow;

            double threshold = 10; // getting expected number presynaptic within threshold is when the approximation stops

            double desiredExpNumPre = desiredNumPrePeak1 * egalPlains[numInputs[post]];

            if (expNumPreHigh < desiredExpNumPre + threshold){
                System.out.println("went high");
                System.out.println(expNumPreHigh);
                System.out.println(desiredExpNumPre);
                System.out.println(egalPlains[numInputs[post]]);

                connectivityWithEgal = CHigh;
            } else if (expNumPreLow > desiredExpNumPre - threshold) {
                connectivityWithEgal = CLow;
            } else {
                double CMid = (CHigh+CLow) / 2.0;
                double expNumPreMid = 0;
                for (int i=1; i<numP; i++){
                    expNumPreMid += min(1.0,CMid * allP[i]);
                }
                while (expNumPreMid > desiredExpNumPre+threshold ||
                       expNumPreMid < desiredExpNumPre-threshold){
                    if (expNumPreMid < desiredExpNumPre){
                        CLow = CMid;
                    } else {
                        CHigh = CMid;
                    }
                    CMid = (CLow+CHigh)/2.0;
                    expNumPreMid = 0;
                    for (int i = 1; i < numP; i++) {
                        expNumPreMid += min(1.0,CMid * allP[i]);
                    }
                }
                connectivityWithEgal = CMid;
                //debugCount++;
                //if (true) System.out.println("expNumPreLow="+expNumPreLow+" desiredExpNumPre="+desiredExpNumPre+" connectivity="+connectivityWithEgal+" expNumPreMid"+expNumPreMid);

            }
        }else{
            connectivityWithEgal = egalPlains[numInputs[post]];//(numInputs[post] < egalHelper.length ? egalHelper[numInputs[post]] + connectivity : connectivity);
        }

        /* actually choose connections now that connectivity is established */
        for (int j = 0; j < sizePreGroup; j++) {
            double sds = preLocs[j].distance(locPost) / lenScale;
            if (sds < 2.0 && post != j && random() < connectivityWithEgal * exp(-sds * sds / 2)) {
                // then make a connection
                cPostFromPre[j][numPostFromPre[j]] = post; //place on connectivity list
                numPostFromPre[j]++; //update number of E's postsynaptic to j
                cPreToPost[post][numPreToPost[post]] = j; //place on presynaptic connectivity list
                numPreToPost[post]++; //update number of E's presynaptic to i
                }
        }
        /*== chose presynaptic connections ==*/
    }

	/* get LGN connectivity data from MATLAB .mat file */
    public final void getLGNInputs(){
        try {
            MatFileReader matfilereader;
            matfilereader = new MatFileReader(LGNInputFile);


            /*==Import out-degrees and out-connectivities
             *of LGN cells (these will need further processing,
             *like changing to java indexing starting at 0)==*/
            double[][] nELGNONPre = ((MLDouble)matfilereader.getMLArray(       "numLGNONToECortex")).getArray();
            double[][] nELGNOFFPre = ((MLDouble)matfilereader.getMLArray(      "numLGNOFFToECortex")).getArray();
            double[][] nILGNONPre = ((MLDouble)matfilereader.getMLArray(       "numLGNONToICortex")).getArray();
            double[][] nILGNOFFPre = ((MLDouble)matfilereader.getMLArray(      "numLGNOFFToICortex")).getArray();
            double[][] cELGNONPre = ((MLDouble)matfilereader.getMLArray(       "LGNONToECortex")).getArray();
            double[][] cELGNOFFPre = ((MLDouble)matfilereader.getMLArray(      "LGNOFFToECortex")).getArray();
            double[][] cILGNONPre = ((MLDouble)matfilereader.getMLArray(       "LGNONToICortex")).getArray();
            double[][] cILGNOFFPre = ((MLDouble)matfilereader.getMLArray(      "LGNOFFToICortex")).getArray();
            double[][] maxNumLGNToCortex = ((MLDouble)matfilereader.getMLArray("maxNumLGNToCortex")).getArray();
            double[][] corticalLocLGNON = ((MLDouble)matfilereader.getMLArray( "corticalLocLGNON")).getArray();
            double[][] corticalLocLGNOFF = ((MLDouble)matfilereader.getMLArray("corticalLocLGNOFF")).getArray();
            double[][] numInputsExcPre = ((MLDouble)matfilereader.getMLArray(  "numInputsExc")).getArray();
            double[][] numInputsInhPre = ((MLDouble)matfilereader.getMLArray(  "numInputsInh")).getArray();
            /*==Imported, now process them==*/
            int maxOutDegree = (int)maxNumLGNToCortex[0][0];
            cELGNON = new int[numLGN][maxOutDegree];
            cELGNOFF = new int[numLGN][maxOutDegree];
            cILGNON = new int[numLGN][maxOutDegree];
            cILGNOFF = new int[numLGN][maxOutDegree];
            nELGNON = new int[numLGN];
            nELGNOFF = new int[numLGN];
            nILGNON = new int[numLGN];
            nILGNOFF = new int[numLGN];

            for (int i=0;i<numLGN;i++){
                nELGNON[i] = (int)nELGNONPre[i][0];
                nELGNOFF[i] = (int)nELGNOFFPre[i][0];
                nILGNON[i] = (int)nILGNONPre[i][0];
                nILGNOFF[i] = (int)nILGNOFFPre[i][0];
                for (int j=0;j<nELGNON[i];j++){
                    cELGNON[i][j] = (int)cELGNONPre[i][j] - 1;
                }
                for (int j=0;j<nELGNOFF[i];j++){
                    cELGNOFF[i][j] = (int)cELGNOFFPre[i][j] - 1;
                }
                for (int j=0;j<nILGNON[i];j++){
                    cILGNON[i][j] = (int)cILGNONPre[i][j] - 1;
                }
                for (int j=0;j<nILGNOFF[i];j++){
                    cILGNOFF[i][j] = (int)cILGNOFFPre[i][j] - 1;
                }
                locLGNON[i] = new MyPoint(corticalLocLGNON[i][0], corticalLocLGNON[i][1]);
                locLGNOFF[i] = new MyPoint(corticalLocLGNOFF[i][0], corticalLocLGNOFF[i][1]);
            }
            for (int i=0;i<numExcitatory;i++){
                numInputsExc[i] = (int)numInputsExcPre[i][0];
            }
            for (int i=0;i<numInhibitory;i++){
                numInputsInh[i] = (int)numInputsInhPre[i][0];
            }
            /*==Processing LGN->cortex connection complete==*/
        } catch (IOException ex) {
            Logger.getLogger(Sim.class.getName()).log(Level.SEVERE, null, ex);
            System.exit(1);
        }
    }

    /* load layer 6 */
    public final void loadFB(){
        try {
            MatFileReader matfilereader;
            matfilereader = new MatFileReader(fbInputFile);
            
            /*==Import out-degrees and out-connectivities
             *of LGN cells (these will need further processing,
             *like changing to java indexing starting at 0)==*/
            double[][] fbSpikeTimesPre = ((MLDouble)matfilereader.getMLArray("fbSpikeTimes")).getArray();
            double[][] fbSpikesPre = ((MLDouble)matfilereader.getMLArray("fbSpikes")).getArray();
            double[][] numFBSpikesPre = ((MLDouble)matfilereader.getMLArray("numFBSpikes")).getArray();
            double[][] fbLocPre = ((MLDouble)matfilereader.getMLArray( "fbLoc")).getArray();
            double[][] numFBPre = ((MLDouble)matfilereader.getMLArray( "numFB")).getArray();

            /* Process */
            numFB = (int)numFBPre[0][0];
            numFBSpikes = (int)numFBSpikesPre[0][0];
            fbSpikeTimes = new double[numFBSpikes+1];
            fbSpikes = new int[numFBSpikes];
            for(int i=0;i<numFBSpikes;i++){fbSpikeTimes[i]=fbSpikeTimesPre[i][0];fbSpikes[i]=(int)fbSpikesPre[i][0];}
            fbSpikeTimes[numFBSpikes] = fbEnd + 1.0;
            locFB = new double[numFB][2];
            for(int i=0;i<numFB;i++){locFB[i][0]=fbLocPre[i][0];locFB[i][1]=fbLocPre[i][1];}

            /* setup the p-value arrays */
            fbInitProb = new double[numFB];
            fbExtraProb = new double[numFB];

            /* allocate fb connectivity matrix */
            numExcPostFB = new int[numFB];
            excPostFB = new int[numFB][400];

            /*==Imported, now create layer6->4C connectivity matrix */
            for (int i=0;i<numExcitatory;i++){
                /* Let user know how much progress is done */
                if (i % 3000 == 0){
                    System.out.println("Processed " + i + " neurons");
                }

                /* find the peak probability from the egalitarian connectivity */
                double peak = .15 * fbEgalMod[numInputsExc[i]];



                /* find presynaptic layer 6 cells */
                double locX = locExc[i].x;
                double locY = locExc[i].y;
                for (int j=0;j<numFB;j++){
                    double dxSqr = sqr(locX - locFB[j][0]);
                    double dySqr = sqr(locY - locFB[j][1]);
                    double dSqr = (dxSqr + dySqr) / (fbESD*fbESD);
                    //System.out.println("neuron (" + locX + ", " + locY + ") and fb ("+locFB[j][0]+", "+locFB[j][1]+") dist sqr normalized is " + dSqr);
                    //try{Thread.sleep(250);}catch(InterruptedException e){System.out.println("caught");}
                    if (dSqr < 2.25){
                        double probPreSyn = peak * exp(-dSqr/2);
                        //System.out.println("prob="+probPreSyn);
                        if (random() < probPreSyn){
                            /* make a connection from jth FB cell to ith E cell*/
                            //System.out.println("success");
                            excPostFB[j][numExcPostFB[j]] = i;
                            numExcPostFB[j]++;
                        }
                    }
                }

                /*
                int numListSoFar = 0;
                int[] listSoFar = new int[0];
                for (int j=0;j<numPreSyn;j++){
                    int nextTrain = preSynList[j];
                    int nextTrainSize = (int)numFBSpikeTimes[nextTrain][0];
                    int nextListSize = numListSoFar + nextTrainSize;
                    int[] nextListSoFar = new int[nextListSize];
                    // merge listSoFar with the next list and place in nextListSoFar
                    int m = 0;
                    int n = 0;
                    int k = 0;
                    while (m < numListSoFar && n < nextTrainSize){
                        if (listSoFar[m] < (int)fbSpikeTimes[nextTrain][n]){
                            nextListSoFar[k] = listSoFar[m];
                            k++;
                            m++;
                        } else {
                            nextListSoFar[k] = (int)fbSpikeTimes[nextTrain][n];
                            k++;
                            n++;
                        }
                    }
                    for (int l=m;l<numListSoFar;l++){
                        nextListSoFar[k] = listSoFar[l];
                        k++;
                    }
                    for (int l=n;l<nextTrainSize;l++){
                        nextListSoFar[k] = (int)fbSpikeTimes[nextTrain][n];
                        k++;
                    }
                    listSoFar = nextListSoFar;
                    numListSoFar = nextListSize;
                }


                for (int j=0;j<numListSoFar;j++){
                    fbSpikeTimesExc[i][j] = listSoFar[j];
                }
                numFBSpikeTimesExc[i] = numListSoFar;
                fbSpikeTimesExc[i][numListSoFar] = fbEnd+1; // To stop firing at end of recording*/
            } //run through E cells

            /* allocate fb connectivity matrix to I's */
            numInhPostFB = new int[numFB];
            inhPostFB = new int[numFB][200];

            /* run through I cells */
            for (int i=0;i<numInhibitory;i++){
                /* Let user know how much progress is done */
                if (i % 3000 == 0){
                    System.out.println("Processed " + i + " neurons");
                }

                /* find the peak probability */
                double peak = fbIPeak;



                /* find presynaptic layer 6 cells */
                double locX = locInh[i].x;
                double locY = locInh[i].y;
                for (int j=0;j<numFB;j++){
                    double dxSqr = sqr(locX - fbLocPre[j][0]);
                    double dySqr = sqr(locY - fbLocPre[j][1]);
                    double dSqr = (dxSqr + dySqr)/ (fbISD*fbISD);
                    if (dSqr < 2.25){
                        double probPreSyn = peak * exp(-dSqr/2);
                        if (random() < probPreSyn){
                            /* make a connection from jth FB cell to ith E cell*/
                            inhPostFB[j][numInhPostFB[j]] = i;
                            numInhPostFB[j]++;
                        }
                    }
                }
            }

            /* allocate space for fbPValues */
            fbPValue = new double[numFB];

            /*==Processing LGN->cortex connection complete==*/
        } catch (IOException ex) {
            Logger.getLogger(Sim.class.getName()).log(Level.SEVERE, null, ex);
            System.exit(1);
        }
    }

    /*==set p-values according to location in cortex==*/
    public void setFBPValues(){
        /* find normalizing factor of convolving Gaussian */
        double SD = 75.0;
        double gridWidth = 30;
        double dx = 2.0 * SD / gridWidth;
        double N = 0.0;
        for (double x = -SD; x < SD; x+=dx){
            for (double y = -SD; y < SD; y+=dx){
                double normDSqr = (sqr(x) + sqr(y)) / sqr(SD);
                if (normDSqr < 1.0) N += exp(-0.5 * normDSqr);
            }
        }
        
        /* go through each layer 6 neuron and assign a p-value */
        for (int i=0;i<numFB;i++){
            /* get the neuron location */
            double[] loc = locFB[i];

            /* do a cheap convolution */
            double conv = 0.0;
            for (double x = loc[0]-SD; x < loc[0]+SD; x+=dx){
                for (double y = loc[1]-SD; y < loc[1]+SD; y+=dx){
                    double dlocx = x - loc[0];
                    double dlocy = y - loc[1];
                    double normDSqr = (sqr(dlocx) + sqr(dlocy)) / sqr(SD);
                    if (normDSqr < 1.0){
                        double[] newLoc = {x, y};
                        conv += constantPValues(newLoc) * exp(-0.5 * normDSqr);
                    }
                }
            }
            conv /= N;

            /* store the result in the current layer 6 neuron */
            fbPValue[i] = conv;
        }

        for (int i=0;i<numFB;i++){
            fbInitProb[i] = fbPValue[i];
            fbExtraProb[i]= fbPValue[i] - 1.0;
        }
    }

    public double constantPValues(double[] loc){
        double constantP;
        /* switch to 1st HC coordinates */
        double[] locHC1 = { NetHelper.movCoordToHC1(loc[0], HCWidth),
                            NetHelper.movCoordToHC1(loc[1], HCWidth)};

        double locLen = MyPoint.distance(0, 0, locHC1[0], locHC1[1]);

        /* find constant value based on angle */
        double angle = NetHelper.pointToAngle(locHC1[0],locHC1[1]);
        //double orientAngle = NetHelper.pointToAngle(kOrient.x,kOrient.y);

        int ind = ((int)floor((angle+2.0*kAngle) / (PI/6.0))) % 12;        

        constantP = pValues[ind];

        return constantP;
    }

	/* set LGN strength */
    public void setSLGN(double newSLGN){
        sLGN = newSLGN;
        sLGNInh = newSLGN;
        sEE = 0.5 * newSLGN;
    }

    public void setSEE(double newSEE){
        sEE = newSEE;
        sLGN = 2 * newSEE;
        sLGNInh = 2 * newSEE;
        sStrongAmbE = newSEE;
    }

    public void setJustSEE(double justSEE){
        this.justSEE = justSEE;
        sStrongAmbE = justSEE;
    }

	/* set the ratio of SEI to SEE */
    public void setSEIOverSEE(double newSEIOverSEE){
        sEI = newSEIOverSEE * sEE;
        sII = newSEIOverSEE * sEE;
    }

    public void setSIEOverSEE(double newSIEOverSEE){
        sIE = newSIEOverSEE * sEE;
        sStrongAmbI = sIE;
    }

	/* set the ambient drive parameters */
    public void setAmbient(double strongRateE, double weakRateE, double strongRateI, double weakRateI){
        rateStrongAmbE = strongRateE;
        rateWeakAmbE = weakRateE;
//        rateTotalE = strongRateE + weakRateE;
//        strongProbE = (rateStrongAmbE+rateWeakAmbE>0 ? rateStrongAmbE / (rateStrongAmbE + rateWeakAmbE) : 1.0);

        rateStrongAmbI = strongRateI;
        rateWeakAmbI = weakRateI;
//        rateTotalI = strongRateI + weakRateI;
//        strongProbI = (rateStrongAmbI+rateWeakAmbI>0 ? rateStrongAmbI / (rateStrongAmbI + rateWeakAmbI) : 1.0);

        ensureFB();
    }

	/* set the I threshold range */
    public void setThreshIRange(double low, double high){
        for (int i=0;i<numInhibitory;i++){
            threshI[i] = low + (high-low)*random();
        }
    }

    public void setThreshERange(double low, double high){
        double mean,std,topCutoff,bottomCutoff;
        mean = (low+high)/2;
        std  = (high-low)/2;
        topCutoff = mean + std;
        bottomCutoff = mean - std;

        for (int i=0;i<numExcitatory;i++){
            threshE[i] = min( topCutoff,
                              max( bottomCutoff,
                                   mean + std * normal.nextDouble()
                                 )
                            );
        }
    }

    public void setSModSeeds(){
        for (int i=0;i<numExcitatory;i++){
            sModSeedExc[i] = random();
        }
        for (int i=0;i<numInhibitory;i++){
            sModSeedInh[i] = random();
        }
        for (int i=0;i<numLGN;i++){
            sModSeedLGNON[i] = random();
            sModSeedLGNOFF[i]= random();
        }
    }

	/* set the range of possible I to I strength modifiers */
    public void setModIIRange(double low, double high){
        modIILow = low;
        modIIHigh = high;
        modIIRange = high - low;
    }

	/* set light grating contrast */
    public void setContrast(double newContrast){
        eps = newContrast;

        ensureFB();
    }

    public void ensureFB(){
        fbAvg = 1.0;// + eps * (fbPref + fbOrth) / 2;
        fbMod = eps;

        setFBPValues();
    }

    public void setIntensity(double newIntensity){
        I_0 = newIntensity;
    }

    public void setKFreq(double newKFreq){
        kFreq = newKFreq;
        setPhases();
        setSpFreqMod();
    }
    public void setPhases(){
        for (int i=0;i<numLGN;i++){
            phaseON[i] = kFreq * kOrient.innerProduct(locLGNON[i]);
            phaseOFF[i]= kFreq * kOrient.innerProduct(locLGNOFF[i]) + PI;
        }
    }

    public double getSpFreqMod(double cpd){
        double sa = 0.066*2.0*1.2;
        double sb = 0.093*2.0*1.2;
        return max(
                    0.0,
                    0.95*(1/0.3359)*(1.0*exp(-sa*sa/4.0/PI*cpd*cpd*4*PI*PI) - 0.74*exp(-sb*sb/4.0/PI*cpd*cpd*4*PI*PI))
                  );

    }

    public void setSpFreqMod(){
        // recall kFreq is in units of radians per micrometer
        double cpd = kFreq * 2000.0 / (2*PI);
        spFreqMod = getSpFreqMod(cpd);
    }

    public void setKOrient(MyPoint newKOrient){
        kOrient = newKOrient;
        kAngle = NetHelper.pointToAngle(kOrient.x, kOrient.y);
        setPhases();
    }

    /* temporal frequency of drift grating */
    public void setFrequency(double newW){
        w = newW;
    }

    /* resets the LGN inputs if a different LGN input file is given--also rebuilds the network */
    public void setLGNInputsIfNew(String LGNInputFile){
        if (!this.LGNInputFile.equals(LGNInputFile)){
            this.LGNInputFile = LGNInputFile;

            System.out.println("New LGN Input File detected. Rebuilding Network...");
            getLGNInputs();
            setupCortexNetwork(); //need to also reset the cortical network bc egal can be different now
            setSModSeeds();
        }
    }

    /* calculate the conductances for the next timestep */
    public void evolveAMPAGABA(double tau, double[][] g, double numConductances){
        for(int i = 0; i < numConductances; i++){
            /* main conductance is g[i][0] */
            g[i][0] = exp(-dt/tau) * (g[i][0] + (dt/tau) * (3 * g[i][1] + (dt/tau) * (3 * g[i][2] + (dt/tau) * g[i][3])));
            g[i][1] = exp(-dt/tau) * (g[i][1] + (dt/tau) * (2 * g[i][2] + (dt/tau) * g[i][3]));
            g[i][2] = exp(-dt/tau) * (g[i][2] + (dt/tau) * g[i][3]);
            g[i][3] *= exp(-dt/tau);
        }
    }

    /* set the conductances after being synapsed on */
    public void synapseAMPAGABA(double tau, double[] g, double spaceK, double S){
        double accum = S*spaceK*exp(-dt/tau) / (6.0 * tau);
        g[3] += accum;
        accum *= dt/tau;
        g[2] += accum;
        accum *= dt/tau;
        g[1] += accum;
        accum *= dt/tau;
        g[0] += accum;
    }

    /* same as above except for NMDA */
    public void evolveNMDA(double[][] g, double numConductances){
        for(int i = 0; i < numConductances; i++){
            /* main conductance is g[i][0] */
            g[i][1] *= exp(-dt/tauNMDA1);
            g[i][2] *= exp(-dt/tauNMDA2);
            g[i][0] = g[i][1] - g[i][2];
        }
    }

    public void synapseNMDA(double[] g, double spaceK, double S){
        double coeff = S*spaceK/(tauNMDA1-tauNMDA2);
        g[1] += coeff;
        g[2] += coeff;
    }

    public final double sqr(double x){ return x*x; }

    public final double gaussian2D(double sigma, double normSqr){
        if (normSqr > 4*sigma*sigma){
            return 0.0;
        } else {
            return 1.0 / (2.0 * PI * sqr(sigma)) * exp(-normSqr / (2.0 * sqr(sigma)));
        }
    }

	/* find the voltage for the next timestep */
    public void evolveVExc(){
        for (int i=0;i<numExcitatory;i++){
            if (sleepExc[i] >= 0){
                vExc[i] = vExc[i] * exp(-dt * gLExc) + dt * ( (gLGNExc[i][0]+gAMPAExc[i][0]+gNMDAExc[i][0]) * (eRev - vExc[i]) + gGABAExc[i][0] * (iRev - vExc[i]) );
            }
        }
    }

    public void evolveVInh(){
        for (int i=0;i<numInhibitory;i++){
            if (sleepInh[i] >= 0){
                vInh[i] = vInh[i] * exp(-dt * gLInh) + dt * ( (gLGNInh[i][0]+gAMPAInh[i][0]+gNMDAInh[i][0]) * (eRev - vInh[i]) + gGABAInh[i][0] * (iRev - vInh[i]) );
            }
        }
    }

	/* move the refractory clocks forward one */
    public void evolveSleep(){
        for (int i=0;i<numExcitatory;i++){
            sleepExc[i]++;
        }
        for (int i=0;i<numInhibitory;i++){
            sleepInh[i]++;
        }
    }

    public void evolveLGN(){
        for (int i=0;i<numLGN;i++){
            /* update the voltage */

            double kickDir = 2*round(random()) - 1;
            vLGNON[i] = vLGNON[i] * exp(-dt * leakLGN) + dt  * I(t,phaseON[i]) + kickDir * shotNoise * (double)poisson.nextInt(dt*noiseRate);

            kickDir = 2*round(random()) - 1;
            vLGNOFF[i] = vLGNOFF[i] * exp(-dt * leakLGN) + dt * I(t,phaseOFF[i]) + kickDir * shotNoise * (double) poisson.nextInt(dt*noiseRate);

            /* check if a spike occured */
            if (vLGNON[i] > threshLGN){
                /* spike occured, send it out */
                for (int j=0;j<nELGNON[i];j++){
                    //synapseAMPAGABA(tauAMPA, gLGNExc[cELGNON[i][j]], 1.0, sLGN * modifier(sModSeedLGNON[i],sModSeedExc[i]));
                    int target = cELGNON[i][j];
                    //delayCylinder[target][(spikeDelayIndex + (int)floor(random() * lgnSpikeDelaySteps)) % spikeDelaySteps] += 1.2 * sLGNOverSEEFactor * modifier(sModSeedLGNON[i],sModSeedExc[target]);
                    synapseAMPAGABA(tauAMPA, gLGNExc[target], 1.0, justSEE * sLGNOverSEEFactor * modifier(sModSeedLGNON[i],sModSeedExc[target]));
                }
                for (int j=0;j<nILGNON[i];j++){
                    synapseAMPAGABA(tauAMPA, gLGNInh[cILGNON[i][j]], 1.0, sLGNInh);
                }

                /* reset voltage, no refractory */
                vLGNON[i] = 0;
            }

            if (vLGNOFF[i] > threshLGN){
                for (int j=0;j<nELGNOFF[i];j++){
                    //synapseAMPAGABA(tauAMPA, gLGNExc[cELGNOFF[i][j]], 1.0, sLGN);
                    int target = cELGNOFF[i][j];
                    //delayCylinder[target][(spikeDelayIndex + (int)floor(random() * lgnSpikeDelaySteps)) % spikeDelaySteps] += 1.2 * sLGNOverSEEFactor * modifier(sModSeedLGNOFF[i],sModSeedExc[target]);
                    synapseAMPAGABA(tauAMPA, gLGNExc[target], 1.0, justSEE * sLGNOverSEEFactor * modifier(sModSeedLGNOFF[i],sModSeedExc[target]));
                }
                for (int j=0;j<nILGNOFF[i];j++){
                    synapseAMPAGABA(tauAMPA, gLGNInh[cILGNOFF[i][j]], 1.0, sLGNInh);
                }

                /* reset voltage, no refractory */
                vLGNOFF[i] = 0;
            }
        } /* for all LGN */
    } /* evolve LGN */

    public double I(double t,double phase){
        double r = 2*PI*w*t + phase;
        double s = sin( r );

        if ( (r/PI+0.5) % 2.00 > 1.0 ){
            return I_0 * (1 + eps * spFreqMod * s);
        } else {
            return I_0 * (1 + eps * spFreqMod * signum(s) * sqrt(abs(s)));
        }
    }

    /* find all conductances for next timestep */
    public void evolveConductances(){
        evolveAMPAGABA(tauAMPA, gAMPAExc, numExcitatory);
        evolveAMPAGABA(tauGABA, gGABAExc, numExcitatory);
        evolveNMDA(gNMDAExc, numExcitatory);

        evolveAMPAGABA(tauAMPA, gAMPAInh, numInhibitory);
        evolveAMPAGABA(tauGABA, gGABAInh, numInhibitory);
        evolveNMDA(gNMDAInh, numInhibitory);

        evolveAMPAGABA(tauAMPA, gLGNExc, numExcitatory);
        evolveAMPAGABA(tauAMPA, gLGNInh, numInhibitory);
    }

    /* implement layer 6 recording feedback */
    public void enactFB(){
        /* check if any spikes happened */
        while (eps > 0.5 && fbClock > fbSpikeTimes[currentFBSpike]){
            /* hit every neuron postsynaptic to fbSpike[currentFBSpike] */
            int firingNeuron = fbSpikes[currentFBSpike];
            boolean initFired = random() < fbInitProb[firingNeuron];
            if (initFired){
                boolean extraFired = random() < fbExtraProb[firingNeuron];
                if (extraFired){
                    for (int k = 0; k < numExcPostFB[firingNeuron]; k++) {
                        /* spike on kth postsynaptic neuron */
                        int target = excPostFB[firingNeuron][k];
                        spikeEE(target, 1.0);                        
                        delayCylinder[target][(int) floor(random() * spikeDelaySteps)] += 1.0;                        
                    }
                    for (int k = 0; k < numInhPostFB[firingNeuron]; k++) {
                        /* spike on kth postsynaptic neuron */
                        int target = inhPostFB[firingNeuron][k];
                        spikeIE(target, 2.0);
                    }
                } else {
                    for (int k = 0; k < numExcPostFB[firingNeuron]; k++) {
                        /* spike on kth postsynaptic neuron */
                        int target = excPostFB[firingNeuron][k];
                        spikeEE(target, 1.0);
                    }
                    for (int k = 0; k < numInhPostFB[firingNeuron]; k++) {
                        /* spike on kth postsynaptic neuron */
                        int target = inhPostFB[firingNeuron][k];
                        spikeIE(target, 1.0);
                    }
                }
            }

            /* update currentFBSpike */
            currentFBSpike++;
        }

        /* update the clock */
        fbClock += dt;
        if (fbClock > fbEnd){
            /* reset the clock */
            fbClock = fbStart;
            currentFBSpike = 0;
        }
    }

    /* implement the ambient drive for this timestep */
    public void enactAmbientDrive(){
        for (int i=0;i<numExcitatory;i++){
            //ambient drive part here
            double finalRateStrongE = rateStrongAmbE * fbEgalMod[numInputsExc[i]];
            double strongProbE = finalRateStrongE / (finalRateStrongE + rateWeakAmbE);
            double numStrongSpikes = (double)poisson.nextInt(dt*(rateWeakAmbE + finalRateStrongE));
            if (numStrongSpikes > 0){
                if (random() < strongProbE){
                    //synapseAMPAGABA(tauAMPA, gLGNExc[i], numStrongSpikes, sStrongAmbE);
                    spikeEE(i, numStrongSpikes);
                    //synapseNMDA(gLGNExc[i], 0.34, sStrongAmbE);
                } else {
                    synapseAMPAGABA(tauAMPA, gLGNExc[i], numStrongSpikes, sWeakAmb);
                }
            }
        }
        for (int i=0;i<numInhibitory;i++){
            //ambient drive part here
            //double numStrongSpikes = (double)poisson.nextInt(dt*rateTotal);
            //if (numStrongSpikes > 0) synapseAMPAGABA(tauAMPA, gLGNInh[i], numStrongSpikes, (random() < strongProb ? sStrongAmb : sWeakAmb));
            double finalRateStrongI = rateStrongAmbI;
            double strongProbI = finalRateStrongI / (finalRateStrongI + rateWeakAmbI);
            double numStrongSpikes = (double)poisson.nextInt(dt*(rateWeakAmbI + finalRateStrongI));
            if (numStrongSpikes > 0){
                if (random() < strongProbI){
                    //synapseAMPAGABA(tauAMPA, gLGNInh[i], numStrongSpikes, sStrongAmbI * 0.66);
                    //synapseNMDA(gNMDAInh[i], 0.34, sStrongAmbI);
                    spikeIE(i,numStrongSpikes);
                } else {
                    synapseAMPAGABA(tauAMPA, gLGNInh[i], numStrongSpikes, sWeakAmb);
                }
            }

        }
    }

    /* find any cortical spikes this timestep and implement them */
    public void enactCorticalSpikes(){
        /* prepare the stats vector */
        numSpikesNowExc = 0;
        numSpikesNowInh = 0;

        delayCylinderStep();

        /* detect E population spikes */
        for(int i = 0; i < numExcitatory; i++){
            if( vExc[i] >= threshE[i] ){
                /* stats collecting */
                eSpikes++;
                spikesNowExc[numSpikesNowExc] = i;
                spikeTimesNowExc[numSpikesNowExc] = t;
                numSpikesNowExc++;

                /* tell the patch as well */
                if(membershipExc[i]==1){
                    patchOneFRExcTemp++;
                } else if(membershipExc[i]==2){
                    patchTwoFRExcTemp++;
                }

                delaySpikeE(i);
            }
        }
        /* among I population */
        for(int i = 0; i < numInhibitory; i++){
            if( vInh[i] >= threshI[i] ){
                /* stats collecting */
                iSpikes++;
                spikesNowInh[numSpikesNowInh] = i;
                spikeTimesNowInh[numSpikesNowInh] = t;
                numSpikesNowInh++;

                /* tell the patch as well */
                if(membershipInh[i]==1){
                    patchOneFRInhTemp++;
                } else if(membershipInh[i]==2){
                    patchTwoFRInhTemp++;
                }

                fastSpikeI(i);
            }
        }
    }

    /* slightly faster version of above */
    public void enactCorticalSpikesNOSTAT(){
        delayCylinderStep();

        /* detect E population spikes */
        for(int i = 0; i < numExcitatory; i++){
            if( vExc[i] >= threshE[i] ){
                /* notify feedback */
                //for(int j=0;j<50;j++) timeLocalFR[j]++; //so timeLocalFR is measured in spikes / 50ms

                delaySpikeE(i);
            }
        }
        /* among I population */
        for(int i = 0; i < numInhibitory; i++){
            if( vInh[i] >= threshI[i] ){
                fastSpikeI(i);
            }
        }
    }

    public void delayCylinderStep(){
        spikeDelayIndex++;
        spikeDelayIndex %= spikeDelaySteps;
        for (int i=0;i<numExcitatory;i++){
            if (delayCylinder[i][spikeDelayIndex] > 0){
                spikeEE(i,delayCylinder[i][spikeDelayIndex]);
                delayCylinder[i][spikeDelayIndex] = 0;
            }
        }
    }

    /* make E cell send delayed spikes */
    public void delaySpikeE(int index){
        /* begin refractory period */
        vExc[index] = 0;
        sleepExc[index] = -(refractoryE / dt);
        /* enact spike */
        /*===For E->E, put on all target cell delay rings===*/
        for (int j = 0; j < numEPostExc[index]; j++) {
            int target = cEExc[index][j];
            delayCylinder[target][(int)floor(random() * spikeDelaySteps)] += modifier(sModSeedExc[index],sModSeedExc[target]);
        }
        for (int j = 0; j < numIPostExc[index]; j++) {
            int target = cIExc[index][j];
            spikeIE(target,1.0);
        }
    }

    /* make I cell send spikes */
    public void fastSpikeI(int index){
        /* begin refractory period */
        vInh[index] = 0;
        sleepInh[index] = -(refractoryI / dt);
        /* enact spike to E's */
        for (int j = 0; j < numEPostInh[index]; j++) {
            int target = cEInh[index][j];
            fastSpikeEI(target,index);
        }

        /* enact spike to I's */
        for (int j = 0; j < numIPostInh[index]; j++) {
            int target = cIInh[index][j];
            fastSpikeII(target);
        }
    }

    public void spikeEE(int target, double modifier){
        synapseAMPAGABA(tauAMPA, gAMPAExc[target], 0.8, justSEE * modifier);
        synapseNMDA(gNMDAExc[target], 0.2, justSEE * modifier);
    }

    public void spikeIE(int target, double modifier){
        synapseAMPAGABA(tauAMPA, gAMPAInh[target], 2.0/3.0, sIE * modifier);
        synapseNMDA(gNMDAInh[target], 1.0/3.0, sIE * modifier);
    }

    public void fastSpikeEI(int target, int source){
        double modifier = sEI * modifier(sModSeedExc[target],sModSeedInh[source]);

        gGABAExc[target][0] += modifier * fastSpikeInitial[3];
        gGABAExc[target][1] += modifier * fastSpikeInitial[2];
        gGABAExc[target][2] += modifier * fastSpikeInitial[1];
        gGABAExc[target][3] += modifier * fastSpikeInitial[0];
    }

    public void fastSpikeII(int target){
        double modifier = (modIILow + modIIRange*random()) * sII;

        gGABAInh[target][0] += modifier * fastSpikeInitial[3];
        gGABAInh[target][1] += modifier * fastSpikeInitial[2];
        gGABAInh[target][2] += modifier * fastSpikeInitial[1];
        gGABAInh[target][3] += modifier * fastSpikeInitial[0];
    }

    public void maintainPatchStats(){
        statStep = (statStep + 1) % statSkipInSteps;
        if (statStep == 0){
            patchOneFRExc = patchOneFRExcTemp;
            patchOneFRExcTemp = 0;

            patchTwoFRExc = patchTwoFRExcTemp;
            patchTwoFRExcTemp = 0;

            patchOneFRInh = patchOneFRInhTemp;
            patchOneFRInhTemp = 0;

            patchTwoFRInh = patchTwoFRInhTemp;
            patchTwoFRInhTemp = 0;
        }
    }

    public void resetStats(){
        eSpikes = 0;
        iSpikes = 0;
    }

    public void resetNetState(){
        resetStats();

        for (int i=0;i<numExcitatory;i++){
            vExc[i] = random() * .9;
            sleepExc[i] = 0;
        }
        for (int i=0;i<numInhibitory;i++){
            vInh[i] = random() * .9;
            sleepInh[i] = 0;
        }
        gLGNExc = new double[numExcitatory][4];
        gAMPAExc = new double[numExcitatory][4];
        gGABAExc = new double[numExcitatory][4];
        gNMDAExc = new double[numExcitatory][3];
        gNMDAInh = new double[numInhibitory][3];
        gGABAInh = new double[numInhibitory][4];
        gLGNInh = new double[numInhibitory][4];
        gAMPAInh = new double[numInhibitory][4];

        /* reset patch stats */
        statStep = 0; //will be incremented mod statStkipsInSteps every step
        patchOneFRExc = 0;
        patchOneFRExcTemp = 0;
        patchOneFRInh = 0;
        patchOneFRInhTemp = 0;
        patchTwoFRExc = 0;
        patchTwoFRExcTemp = 0;
        patchTwoFRInh = 0;
        patchTwoFRInhTemp = 0;

        /* reset e-delay cylinder */
        for (int i=0;i<numExcitatory;i++){
            for (int j=0;j<spikeDelaySteps;j++) delayCylinder[i][j] = 0.0;
        }

        t = 0;
        numSteps=0;
    }

    /*===Go forward in time one time step===*/
    public void step(){
        evolveVExc();
        evolveVInh();
        evolveConductances();
        evolveSleep();
        evolveLGN();

        enactAmbientDrive();
        enactFB();
        enactCorticalSpikes();

        maintainPatchStats();

        /* update time */
        t += dt;
    }

    /*===Step forward without recording===*/
    public void stepNOSTAT(){
        evolveVExc();
        evolveVInh();
        evolveConductances();
        evolveSleep();
        evolveLGN();

        enactAmbientDrive();
        enactFB();
        enactCorticalSpikesNOSTAT();

        /* update time */
        t += dt;
    }

    public void jump(int numTimeSteps){
        for (int i=0;i<numTimeSteps;i++){
            step();
        }
    }

    public void jumpNOSTAT(int numTimeSteps){
        for (int i=0;i<numTimeSteps;i++){
            stepNOSTAT();
        }
    }

    /* the following is for saving and loading data to disk */
    public void saveFBConnectivity(String fname){
        try{
            FileOutputStream fos = new FileOutputStream(fname);
            ObjectOutputStream oos = new ObjectOutputStream(fos);
            oos.writeInt(numFB);
            oos.writeObject(excPostFB);
            oos.writeObject(inhPostFB);
            oos.writeObject(numExcPostFB);
            oos.writeObject(numInhPostFB);
            oos.writeObject(locFB);
        } catch(Exception e){}
    }

    public void loadFBConnectivity(String fname){
        try{
            FileInputStream fis = new FileInputStream(fname);
            ObjectInputStream ois = new ObjectInputStream(fis);
            numFB = ois.readInt();
            excPostFB = (int[][]) ois.readObject();
            inhPostFB = (int[][]) ois.readObject();
            numExcPostFB = (int[]) ois.readObject();
            numInhPostFB = (int[]) ois.readObject();
            locFB = (double[][]) ois.readObject();
        } catch(Exception e){}
    }

    public void saveConnectivity(String fname){
        try{
            FileOutputStream fos = new FileOutputStream(fname);
            ObjectOutputStream oos = new ObjectOutputStream(fos);
            oos.writeObject(numEPostExc);
            oos.writeObject(numIPostExc);
            oos.writeObject(numEPostInh);
            oos.writeObject(numIPostInh);
            oos.writeObject(numEPreExc);
            oos.writeObject(numIPreExc);
            oos.writeObject(numEPreInh);
            oos.writeObject(numIPreInh);
            oos.writeObject(cEExc);
            oos.writeObject(cIExc);
            oos.writeObject(cEInh);
            oos.writeObject(cIInh);
            oos.writeObject(cPreEExc);
            oos.writeObject(cPreIExc);
            oos.writeObject(cPreEInh);
            oos.writeObject(cPreIInh);

        } catch(Exception e){}
    }

    public void loadConnectivity(String fname){
        try{
            FileInputStream fis = new FileInputStream(fname);
            ObjectInputStream ois = new ObjectInputStream(fis);
            numEPostExc = (int[])ois.readObject();
            numIPostExc = (int[])ois.readObject();
            numEPostInh = (int[])ois.readObject();
            numIPostInh = (int[])ois.readObject();
            numEPreExc = (int[])ois.readObject();
            numIPreExc = (int[])ois.readObject();
            numEPreInh = (int[])ois.readObject();
            numIPreInh = (int[])ois.readObject();
            cEExc = (int[][])ois.readObject();
            cIExc = (int[][])ois.readObject();
            cEInh = (int[][])ois.readObject();
            cIInh = (int[][])ois.readObject();
            cPreEExc = (int[][])ois.readObject();
            cPreIExc = (int[][])ois.readObject();
            cPreEInh = (int[][])ois.readObject();
            cPreIInh = (int[][])ois.readObject();
        } catch(Exception e){}
    }
}