const int MEMORYDEPTH = 2500;     //Number of samples per waveform
const int N_WF = 10000;             //Number of WF in file
const int RES = 14;                //Resolution of the digitizer
const int PREPULSE_TICKS = 935;   //Template pre-pulse ticks
const double TICK_LEN = 0.004;    //In mu_s

const double SAT_UP = 19.;
const double SAT_LOW= -19.;

//Calibration
const int INT_LOW = 975;         //Lower limit of wf integration
const int INT_UP  = 975+450;         //Upper fino a 3000
const int NBINS     = 2000;       // " " number of bins
const int NMAXPEAKS = 6;          // " " number of expected peaks
const int MU0_LOW = -300.;     //
const int MU0_UP = 400.;       //This is the cut under which an event is classified as noise
const int SPE_LOW = 1200;  // 0: 5000->7800
const int SPE_UP = 1800;   // 1: 4200->6800
const int S0_LOW = 30;      //
const int S0_UP = 500;
const int SC_LOW = 20;      //
const int SC_UP = 380;
const double FIT_LOW = -1000;
const double FIT_UP = 15000;
const double HMIN = -2000;
const double HMAX = 12000;
//Deconvolution
const int INT_PROMPT = 7500;
const double F_PROMPT = 0.8;
/*
 const double C_FR = 1.;
 const double SPE_AMPL = 7.1;
 */

std::string WF_FILE ("C3_ch1_1953.dat");

