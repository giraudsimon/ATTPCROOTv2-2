#define cRED "\033[1;31m"
#define cYELLOW "\033[1;33m"
#define cNORMAL "\033[0m"
#define cGREEN "\033[1;32m"

bool check_file(const std::string& name);

void run_unpack2(TString dataFile = "runfiles/ar46_run_0085.txt",TString parameterFile = "ATTPC.e15503b.par",TString mappath="/data/ar46/run_0085/")
{

  if(!check_file(dataFile.Data())){
    std::cout<<cRED<<" Run file "<<dataFile.Data()<<" not found! Terminating..."<<cNORMAL<<std::endl;
    exit(0);
  }



  // -----   Timer   --------------------------------------------------------
 TStopwatch timer;
 timer.Start();
 // ------------------------------------------------------------------------

  gSystem->Load("libXMLParser.so");
  // -----------------------------------------------------------------
  // Set file names
  TString scriptfile = "Lookup20150611.xml";
  TString dir = getenv("VMCWORKDIR");
  TString scriptdir = dir + "/scripts/"+ scriptfile;
  TString dataDir = dir + "/macro/data/";
  TString geomDir = dir + "/geometry/";
  gSystem -> Setenv("GEOMPATH", geomDir.Data());

  //TString inputFile   = dataDir + name + ".digi.root";
  //TString outputFile  = dataDir + "output.root";
  TString outputFile  = "output.root";
  //TString mcParFile   = dataDir + name + ".params.root";
  TString loggerFile  = dataDir + "ATTPCLog.log";
  TString digiParFile = dir + "/parameters/" + parameterFile;
  TString geoManFile  = dir + "/geometry/ATTPC_v1.2.root";

  TString inimap   = mappath + "inhib.txt";
  TString lowgmap  = mappath + "lowgain.txt";
  TString xtalkmap = mappath + "beampads_e15503b.txt";

  // -----------------------------------------------------------------
  // Logger
  FairLogger *fLogger = FairLogger::GetLogger();
  fLogger -> SetLogFileName(loggerFile);
  fLogger -> SetLogToScreen(kTRUE);
  fLogger -> SetLogToFile(kTRUE);
  fLogger -> SetLogVerbosityLevel("LOW");

  FairRunAna* run = new FairRunAna();
  run -> SetOutputFile(outputFile);
  //run -> SetGeomFile("../geometry/ATTPC_Proto_v1.0.root");
  run -> SetGeomFile(geoManFile);

  FairRuntimeDb* rtdb = run->GetRuntimeDb();
  FairParAsciiFileIo* parIo1 = new FairParAsciiFileIo();
  parIo1 -> open(digiParFile.Data(), "in");
  //FairParRootFileIo* parIo2 = new FairParRootFileIo();
  //parIo2 -> open("param.dummy_proto.root");
 // rtdb -> setFirstInput(parIo2);
  rtdb -> setSecondInput(parIo1);

  // Settings
  Bool_t fUseDecoder = kTRUE;
  if (dataFile.IsNull() == kTRUE)
    fUseDecoder = kFALSE;

  Bool_t fUseSeparatedData = kFALSE;
    if (dataFile.EndsWith(".txt"))
      fUseSeparatedData = kTRUE;

 /*
  *     Unpacking options:
  *         - SetUseSeparatedData:      To be used with 10 CoBo files without merging. Mainly for the ATTPC. Enabled if the input file is a txt.
  *         - SetPseudoTopologyFrame:   Used to force the graw file to have a Topology frame.
  *         - SetPersistance:           Save the unpacked data into the root file.
  *         - SetMap:                   Chose the lookup table.
  *         - SetMapOpt                 Chose the pad plane geometry. In addition forces the unpacker to use Basic Frames for 1 single file (p-ATTPC case) of Layered
  *                                     Frames for Merged Data (10 Cobos merged data).
  */
  ATDecoder2Task *fDecoderTask = new ATDecoder2Task();
  fDecoderTask -> SetUseSeparatedData(fUseSeparatedData);
  if(fUseSeparatedData) fDecoderTask -> SetPseudoTopologyFrame(kTRUE);//! This calls the method 10 times so for less than 10 CoBos ATCore2 must be modified
  //fDecoderTask -> SetPositivePolarity(kTRUE);
  fDecoderTask -> SetPersistence(kFALSE);
  fDecoderTask -> SetMap(scriptdir.Data());
  fDecoderTask -> SetNumCobo(10);
  fDecoderTask -> SetInhibitMaps(inimap,lowgmap,xtalkmap); // TODO: Only implemented for fUseSeparatedData!!!!!!!!!!!!!!!!!!!1
  fDecoderTask -> SetMapOpt(0); // ATTPC : 0  - Prototype: 1 |||| Default value = 0


  if (!fUseSeparatedData)
    fDecoderTask -> AddData(dataFile);
  else {
    std::ifstream listFile(dataFile.Data());
    TString dataFileWithPath;
    Int_t iCobo = 0;
    while (dataFileWithPath.ReadLine(listFile)) {

      if(!check_file(dataFileWithPath.Data())){
        std::cout<<cRED<<" GRAW file "<<dataFileWithPath.Data()<<" not found! Terminating..."<<cNORMAL<<std::endl;
        exit(0);
      }else{

              if (dataFileWithPath.Contains(Form("CoBo%i",iCobo)) )
                    fDecoderTask -> AddData(dataFileWithPath, iCobo);
              else{
                iCobo++;
                fDecoderTask -> AddData(dataFileWithPath, iCobo);
              }
      }

    }
  }

  run -> AddTask(fDecoderTask);

  ATPSATask *psaTask = new ATPSATask();
  psaTask -> SetPersistence(kTRUE);
  psaTask -> SetThreshold(20);
  psaTask -> SetPSAMode(1); //NB: 1 is ATTPC - 2 is pATTPC
	//psaTask -> SetPeakFinder(); //NB: Use either peak finder of maximum finder but not both at the same time
	psaTask -> SetMaxFinder();
  psaTask -> SetBaseCorrection(kTRUE); //Directly apply the base line correction to the pulse amplitude to correct for the mesh induction. If false the correction is just saved
  psaTask -> SetTimeCorrection(kTRUE); //Interpolation around the maximum of the signal peak. Only affect Z calibration at PSA stage
  run -> AddTask(psaTask);

  ATHoughTask *HoughTask = new ATHoughTask();
	HoughTask ->SetPersistence();
	//HoughTask ->SetLinearHough();
	HoughTask ->SetCircularHough();
  HoughTask ->SetHoughThreshold(100.0); // Charge threshold for Hough
  HoughTask ->SetEnableMap(); //Enables an instance of the ATTPC map:  This enables the MC with Q instead of position
  HoughTask ->SetMap(scriptdir.Data());
	run -> AddTask(HoughTask);


  run -> Init();

  //run -> RunOnTBData();
  run->Run(0,50);

  std::cout << std::endl << std::endl;
  std::cout << "Macro finished succesfully."  << std::endl << std::endl;
  std::cout << "- Output file : " << outputFile << std::endl << std::endl;
  // -----   Finish   -------------------------------------------------------
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << endl;
  cout << endl;
  // ------------------------------------------------------------------------

  gApplication->Terminate();

}

bool check_file(const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

/*fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo0_run_0122_14-08-15_13h27m28s.graw",0);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo0_run_0122_14-08-15_13h27m28s.1.graw",0);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo0_run_0122_14-08-15_13h27m28s.2.graw",0);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo0_run_0122_14-08-15_13h27m28s.3.graw",0);

fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo1_run_0122_14-08-15_13h27m28s.graw",1);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo1_run_0122_14-08-15_13h27m28s.1.graw",1);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo1_run_0122_14-08-15_13h27m28s.2.graw",1);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo1_run_0122_14-08-15_13h27m28s.3.graw",1);

fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo2_run_0122_14-08-15_13h27m29s.graw",2);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo2_run_0122_14-08-15_13h27m29s.1.graw",2);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo2_run_0122_14-08-15_13h27m29s.2.graw",2);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo2_run_0122_14-08-15_13h27m29s.3.graw",2);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo2_run_0122_14-08-15_13h27m29s.4.graw",2);

fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo3_run_0122_14-08-15_13h27m29s.graw",3);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo3_run_0122_14-08-15_13h27m29s.1.graw",3);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo3_run_0122_14-08-15_13h27m29s.2.graw",3);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo3_run_0122_14-08-15_13h27m29s.3.graw",3);

fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo4_run_0122_14-08-15_13h27m29s.graw",4);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo4_run_0122_14-08-15_13h27m29s.1.graw",4);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo4_run_0122_14-08-15_13h27m29s.2.graw",4);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo4_run_0122_14-08-15_13h27m29s.3.graw",4);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo4_run_0122_14-08-15_13h27m29s.4.graw",4);

fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo5_run_0122_14-08-15_13h27m29s.graw",5);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo5_run_0122_14-08-15_13h27m29s.1.graw",5);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo5_run_0122_14-08-15_13h27m29s.2.graw",5);

fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo6_run_0122_14-08-15_13h27m29s.graw",6);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo6_run_0122_14-08-15_13h27m29s.1.graw",6);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo6_run_0122_14-08-15_13h27m29s.2.graw",6);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo6_run_0122_14-08-15_13h27m29s.3.graw",6);

fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo7_run_0122_14-08-15_13h27m29s.graw",7);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo7_run_0122_14-08-15_13h27m29s.1.graw",7);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo7_run_0122_14-08-15_13h27m29s.2.graw",7);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo7_run_0122_14-08-15_13h27m29s.3.graw",7);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo7_run_0122_14-08-15_13h27m29s.4.graw",7);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo7_run_0122_14-08-15_13h27m29s.5.graw",7);

fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo8_run_0122_14-08-15_13h27m29s.graw",8);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo8_run_0122_14-08-15_13h27m29s.1.graw",8);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo8_run_0122_14-08-15_13h27m29s.2.graw",8);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo8_run_0122_14-08-15_13h27m29s.3.graw",8);

fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo9_run_0122_14-08-15_13h27m29s.graw",9);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo9_run_0122_14-08-15_13h27m29s.1.graw",9);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo9_run_0122_14-08-15_13h27m29s.2.graw",9);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo9_run_0122_14-08-15_13h27m29s.3.graw",9);
fDecoderTask -> AddData("/home/ayyadlim/Desktop/ATTPC/run_0122/CoBo9_run_0122_14-08-15_13h27m29s.4.graw",9);*/
