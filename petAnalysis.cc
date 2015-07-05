#include <algorithm> 
#include "barycenterAlgorithm.hh"
#include "TF1.h"
#include<petAnalysis.h>

#define CUT 0

ClassImp(petAnalysis)

//==========================================================================
petAnalysis::petAnalysis(gate::VLEVEL vl, std::string label) : 
IAlgo(vl,"petAnalysis",0,label){
//==========================================================================


}

//==========================================================================
petAnalysis::petAnalysis(const gate::ParamStore& gs, 
			   gate::VLEVEL vl, std::string label) :
  IAlgo(gs,vl,"petAnalysis",0,label){
//==========================================================================


}

//==========================================================================
bool petAnalysis::initialize(){
//==========================================================================

  _m.message("Intializing algorithm",this->getAlgoLabel(),gate::NORMAL);

  return true;  
}

//==========================================================================
bool petAnalysis::execute(gate::Event& evt){
//==========================================================================

  _m.message("Executing algorithm",this->getAlgoLabel(),gate::VERBOSE);
  
  _m.message("Event number:",evt.GetEventID(),gate::VERBOSE);

  // Search primary particle and its first daughter
  gate::MCParticle primary;
  gate::MCParticle firstDaughter;
  findFirstParticle(evt.GetMCParticles(),primary);
  findFirstParticle(primary.GetDaughters(),firstDaughter);

  //True Vertex
  gate::Point3D trueVertex = firstDaughter.GetInitialVtx(); 

 //Try only events with photoelectric and one vertex
  if(firstDaughter.GetCreatorProc() == std::string("phot") 
		  && firstDaughter.GetDaughters().size()==0){

   //   std::cout << "Event number:" << evt.GetEventID() << "\t(" << "x = " << trueVertex.x() << "\ty = "<< trueVertex.y() << "\t z = " << trueVertex.z() << ")" << std::endl; 

	  //Classify sensor hits per planes
	  std::vector<std::vector<gate::Hit*> > planes(6);
	  splitHitsPerPlane(evt,planes);

	  //Apply cut per plane
	  std::vector<std::vector<gate::Hit*> > planesCut(6);
	  for(unsigned int i=0; i<6;i++){
		  applyCut(planes[i], CUT ,planesCut[i]);
	  }

	  std::cout << evt.GetEventID();
	  printSensors(planesCut);
	  std::cout << "\t" << trueVertex.x() << "\t" << trueVertex.y() << "\t" << trueVertex.z() << std::endl;
  }
  return true;
}

//==========================================================================
bool petAnalysis::finalize(){
	//==========================================================================

	_m.message("Finalising algorithm",this->getAlgoLabel(),gate::NORMAL);

  gate::Run* run = &gate::Centella::instance()->getRun();
  int nevt = gate::int_from_string(run->fetch_sstore("num_events"));
  _m.message("Number of generated events in file:",nevt,gate::NORMAL);

  return true;

}

void petAnalysis::findFirstParticle(const std::vector<const gate::MCParticle*> particles, gate::MCParticle& first){
	double firstTime = 10000000.;
	for(unsigned int i=0; i<particles.size();i++){
		//Search the first daughter
		if(particles[i]->GetInitialVtx4D().GetT() < firstTime){
			first = *particles[i];
			firstTime = particles[i]->GetInitialVtx4D().GetT();
		}
	}
}

void petAnalysis::findFirstParticle(std::vector<gate::MCParticle*> particles, gate::MCParticle& first){
	double firstTime = 10000000.;
	for(unsigned int i=0; i<particles.size();i++){
		//Search the first daughter
		if(particles[i]->GetInitialVtx4D().GetT() < firstTime){
			first = *particles[i];
			firstTime = particles[i]->GetInitialVtx4D().GetT();
		}
	}
}

bool petAnalysis::timeOrderParticles(const gate::MCParticle* p1, const gate::MCParticle* p2){
	return (p1->GetInitialVtx4D().GetT() < p2->GetInitialVtx4D().GetT()) ;
}

bool petAnalysis::chargeOrderSensorsDesc(const gate::Hit* s1, const gate::Hit* s2){
	return (s1->GetAmplitude() > s2->GetAmplitude()) ;
}

bool petAnalysis::chargeOrderSensorsAsc(const gate::Hit* s1, const gate::Hit* s2){
	return (s1->GetAmplitude() < s2->GetAmplitude()) ;
}

double petAnalysis::distance(gate::Point3D& p1, gate::Point3D& p2){
	double x = (p1.x() - p2.x());
	double y = (p1.y() - p2.y());
	double z = (p1.z() - p2.z());
	return std::sqrt(std::pow(x,2) + std::pow(y,2) + std::pow(z,2));
}

void petAnalysis::splitHitsPerPlane(gate::Event& evt, std::vector<std::vector<gate::Hit*> >& planes){
	for(unsigned int i=0;i<evt.GetMCSensHits().size(); i++){
		int id = evt.GetMCSensHits()[i]->GetSensorID();
		if(id < 100){
			planes[0].push_back(evt.GetMCSensHits()[i]);
		}else if(id < 2000){
			planes[1].push_back(evt.GetMCSensHits()[i]);
		}else if(id < 3000){
			planes[2].push_back(evt.GetMCSensHits()[i]);
		}else if(id < 4000){
			planes[3].push_back(evt.GetMCSensHits()[i]);
		}else if(id < 5000){
			planes[4].push_back(evt.GetMCSensHits()[i]);
		}else if(id < 6000){
			planes[5].push_back(evt.GetMCSensHits()[i]);
		}
	}
}

// Recieve sensor hits and a percentage of the max count.
// Returns hits with amplitude >= cut*maxAmplitude in filtered
void petAnalysis::applyCut(const std::vector<gate::Hit*>& sensorHits, double cut, std::vector<gate::Hit*>& filtered){
	if(sensorHits.size() > 0){
		gate::Hit* max = *std::max_element(sensorHits.begin(),sensorHits.end(),chargeOrderSensorsAsc);

		for(unsigned int i=0; i<sensorHits.size(); i++){
			if(sensorHits[i]->GetAmplitude() >= cut*max->GetAmplitude()){
				filtered.push_back(sensorHits[i]);
			}
		}
	}
}

void petAnalysis::printSensors(std::vector<std::vector<gate::Hit*> >& planes){
	int id;
	double count=0.;

	for(int i=0;i<100;i++){
		id = i;
		count = findSensors(planes[0],id);
	//	std::cout << "\tid" << id << ": " << count;
		std::cout << "\t" << count;
	}
	for(int i=0;i<100;i++){
		id = 1000 + i;
		count = findSensors(planes[1],id);
	//	std::cout << "\tid" << id << ": " << count;
		std::cout << "\t" << count;
	}
	for(int i=0;i<100;i++){
		id = 2000 + i;
		count = findSensors(planes[2],id);
	//	std::cout << "\tid" << id << ": " << count;
		std::cout << "\t" << count;
	}
	for(int i=0;i<100;i++){
		id = 3000 + i;
		count = findSensors(planes[3],id);
	//	std::cout << "\tid" << id << ": " << count;
		std::cout << "\t" << count;
	}
	for(int i=0;i<100;i++){
		id = 4000 + i;
		count = findSensors(planes[4],id);
	//	std::cout << "\tid" << id << ": " << count;
		std::cout << "\t" << count;
	}
	for(int i=0;i<100;i++){
		id = 5000 + i;
		count = findSensors(planes[5],id);
	//	std::cout << "\tid" << id << ": " << count;
		std::cout << "\t" << count;
	}
}

double petAnalysis::findSensors(std::vector<gate::Hit*>& plane, int id){
	double amplitude = 0;
	for(unsigned int i=0;i<plane.size();i++){
		if(plane[i]->GetSensorID() == id){
			amplitude = plane[i]->GetAmplitude();
	//		std::cout << "(" << plane[i]->GetPosition().x() << "," << plane[i]->GetPosition().y() << "," << plane[i]->GetPosition().z() << ")";
			break;
		}
	}
	return amplitude;
}

bool petAnalysis::chargeOrderPlanesDesc(std::pair<int,double> s1, std::pair<int,double> s2){
	return (s1.second > s2.second);
}

void petAnalysis::computeBarycenters(std::vector<std::vector<gate::Hit*> > planes, std::vector<std::vector<double> >& points, std::vector<std::vector<double> >& errors){
	std::string planesDirections[6] = {"xy","yz","xy","yz","xz","xz"};
	util::barycenterAlgorithm* barycenter = new util::barycenterAlgorithm();

	for(unsigned int i=0;i<6;i++){
		barycenter->setPlane(planesDirections[i]);
		barycenter->computePosition(planes[i]);
		points[i][0] = barycenter->getX1();
		points[i][1] = barycenter->getX2();
		errors[i][0] = barycenter->getX1Err();
		errors[i][1] = barycenter->getX2Err();
	}

}
