////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 
// WriteEvent class is used to create particle list so as to write ROOT and Oscar output 
// WriteTrack stores the information at track level
//  
// Event level functions:
//     	GetNStable / SetNStable()    	- returns/sets the total number of stable particles in an event
//	AddEntry() 		     	- adds a track to the event
//	GetNEntries()			- returns total number of tracks in an event (all of parents and daughters)
//	GetWriteTrack(int i)		- returns the ith track
//
//
// Track level functions:
//	GetFinal()			- returns 1 if particle is stable else 0 if it decays
//	GetNum()			- returns the index of the track
//	GetID()				- returns the PDG ID of the particle
//	GetIst()			- always 0 (required in OSCAR1999A format)
//	GetPx()				- returns x-component of track momentum
//	GetPy()				- returns y-component of track momentum
//	GetPz()				- returns z-component of track momentum
//	GetEnergy()			- returns energy of the track
//	GetMass()			- returns mass of the track
//	GetBranch()			- returns the decay branch where it comes from
//	GetParentIndex()		- returns the index of the immediate parent particle
//	GetWeight()			- returns the weight of the track due to its branching ratio
//		
//	Every Get function comes with a similar Set function in order to set the variables
//
//	OpenOutputFile 			- opens a text file to write out the OSCAR output
//	WriteOscarEvent(std::ofstream* output_file, WriteEvent Event) - writes out OSCAR output directly [refer to WriteOscarOutput.C for more info]
//
//  In order to get a ROOT output refer to WriteROOTOutput.C
//
//
// Roli Esha 05/13/2022
//

#ifndef WriteEvent_h
#define WriteEvent_h

class WriteTrack : public TObject{

   	private:

   	int 	isFinal;
   	int 	num;
   	int 	id;
   	int 	ist;
      	float   px;           
      	float   py;           
      	float   pz;           
      	float   en;
      	float   mass;
      	float   xpos;
      	float   ypos;
      	float 	zpos;
      	double 	br;
      	int 	ParentIndex;
      	double 	weight;

   	public:

      	WriteTrack(){  

      		isFinal = -999;
      		num = -999;
      		id = -999;
      		ist = -999;
	        px = -999;
	        py = -999;
	        pz = -999;
	        en = -999;
	        mass = -999;
	        xpos = -999;
	        ypos = -999;
	        zpos = -999;
	        br  = -999;
	        ParentIndex = -999;
	        weight = -999;
      	};

    	//virtual ~WriteTrack();

    	int  	 GetFinal() {return isFinal; };
		int    	 GetNum() { return num; };
		int      GetID() { return id; };
		int      GetIst() { return ist; };
		float    GetPx() { return px; };
		float    GetPy() { return py; };
		float    GetPz() { return pz; };
		float    GetEnergy() { return en; };
		float    GetMass() { return mass; };
		float    GetXpos() { return xpos; };
		float    GetYpos() { return ypos; };
		float    GetZpos() { return zpos; };
		double      GetBranch() { return br; };
		int      GetParentIndex() { return ParentIndex; };
		double      GetWeight() { return weight; };
	    
        void     SetFinal(int sisFinal) { isFinal = sisFinal;};
        void   	 SetNum(int snum)  	{ num = snum; };
        void     SetID(int sid)  	{ id = sid; };
        void     SetIst(int sist)  	{ ist = sist; };
        void     SetPx(float spx)  	{ px = spx; };
        void     SetPy(float spy)  	{ py = spy; };
        void     SetPz(float spz)  	{ pz = spz; };
        void     SetEnergy(float sen)  { en = sen; };
        void     SetMass(float smass)  { mass = smass; };
        void     SetXpos(float sxpos)  { xpos = sxpos; };
        void     SetYpos(float sypos)  { ypos = sypos; };
        void     SetZpos(float szpos)  { zpos = szpos; };
        void     SetBranch(double sbr)  { br = sbr; };
        void     SetParentIndex(int sParentIndex)  	{ ParentIndex = sParentIndex; };
        void     SetWeight(double sweight)  	{ weight = sweight; };
    
    ClassDef(WriteTrack,1)  
};


class WriteEvent : public TObject {

   	private:

   		int 	stable_particles;
 
   		std::vector<WriteTrack> WriteTrackList;

   	public:	
        WriteEvent(){ stable_particles = 0;}
   		//virtual ~WriteEvent();

      	void    SetNStable(int sstable_particles) { stable_particles = sstable_particles; };
      	int     GetNStable() { return stable_particles; }
      			 	
   		void       ClearEvent();
   		
   		void       AddEntry(WriteTrack newEntry) { WriteTrackList.push_back(newEntry); };

   		Long64_t    GetNEntries() { return WriteTrackList.size(); };
   		WriteTrack  GetWriteTrack(int i) {return WriteTrackList[i]; };

   	ClassDef(WriteEvent,1)  
};

std::ofstream*		OpenOutputFile(const char * file);
void 				WriteOscarEvent(std::ofstream* output_file, WriteEvent Event);

#endif