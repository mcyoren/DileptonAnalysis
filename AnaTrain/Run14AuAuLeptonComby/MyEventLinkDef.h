#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace MyDileptonAnalysis;
#pragma link C++ nestedclass;
#pragma link C++ nestedfunction;
#pragma link C++ nestedtypedef;

#pragma link C++ class MyDileptonAnalysis::MyEventContainer+;
#pragma link C++ class MyDileptonAnalysis::MyEvent+;
#pragma link C++ class MyDileptonAnalysis::MyTrack+;
#pragma link C++ class MyDileptonAnalysis::MyHadron+;
#pragma link C++ class MyDileptonAnalysis::MyElectron+;
#pragma link C++ class MyDileptonAnalysis::MyVTXHit+;
#pragma link C++ class MyDileptonAnalysis::MyPair+;

#pragma link C++ class std::vector<MyDileptonAnalysis::MyTrack>;
#pragma link C++ class std::vector<MyDileptonAnalysis::MyHadron>;
#pragma link C++ class std::vector<MyDileptonAnalysis::MyElectron>;
#pragma link C++ class std::vector<MyDileptonAnalysis::MyVTXHit>;
#pragma link C++ class std::vector<MyDileptonAnalysis::MyEvent>;
#pragma link C++ class std::vector<MyDileptonAnalysis::MyPair>;

#endif /* __CINT__ */
