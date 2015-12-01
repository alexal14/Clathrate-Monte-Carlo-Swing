/*****************************************************************************

    Monte Carlo Simulation of sH Clathrate

    Copyright 20014, 2015 Alexander A. Atamas
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*****************************************************************************/

package model;

import java.awt.Point;
import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Random;

public abstract class MetropolisEnsemble implements IMolecularSimulation {
	
    final double alw_SPCE = 109.47;
    final double dHOw_SPCE = 1.0000;
    
    final double alw_TIP3P = 104.52;
    final double dHOw_TIP3P = 0.9572;
    
    final double alw = alw_TIP3P;
    final double dHOw = dHOw_TIP3P;
    
	protected final int numbMolsTotal = 320;
	protected final int numbWaterMols = 272;
	protected final int numbSmallCages = 24;
	protected final int numbMediumCages = 16;
	protected final int numbHydrogenMols = numbSmallCages + numbMediumCages;
	protected final int numbLargeCages = 8;
	
	protected int stepMC;
	protected boolean go = true;
	protected int freqOfSysOut = 2000;
	
    protected int numbOfMolsActual;
	protected AtomVector box = new AtomVector(24.6717302871656, 21.1679306549331, 19.7305382018905);
	protected double rcutSq;
	protected int numbMCSimulations;	
	
	protected ArrayList<Molecule> molecules = new ArrayList<Molecule>(numbMolsTotal);	
	
	protected ArrayList<Molecule> moleculesH2;
	protected ArrayList<Molecule> moleculesACET;
	protected ArrayList<Molecule> moleculesTHF;	
	protected HashMap<Integer, MolecularModel> molModelMap;
	
	protected double temperature = 300;
	protected double invTemperature = 1/temperature;
	protected double pressure = 1.0;
	protected double volume;
	
	protected double enTot;
	
	protected final double dR = 0.01;
	protected final int ng = 3000; 
	protected double [][] tableEnergyAllMoleciles;
	protected int [] g = new int [ng];
	
	ArrayList<Point2D> gr = new ArrayList<Point2D>();
	
	protected SimulOutputData simulOutputData = new SimulOutputData();
	
	protected MetropolisEnsemble(){
		
	}
	
	protected MetropolisEnsemble(int numbMCSimulations){

		this.numbMCSimulations = numbMCSimulations;		
		initData();
	}
	
	protected MetropolisEnsemble(int numbMCSimulations, double temperature){

		this.numbMCSimulations = numbMCSimulations;	
		this.temperature = temperature;
		initData();
	}
	
	protected MetropolisEnsemble(int numbMCSimulations, double temperature, double pressure){

		this.numbMCSimulations = numbMCSimulations;	
		this.temperature = temperature;
		this.pressure = pressure;
		initData();
	}
	
	protected abstract void engineMonteCarlo(int i, int iMol, int iRan);
	
	private void initData(){
		
		MolecularModel waterModelWater = setupWaterModelTIP3P(alw, dHOw);
//		System.out.println(waterModelWater);
		
		MolecularModel hydrogenModelProm_SM_H2 = setupHydrogenModel();
//		System.out.println(hydrogenModelProm_SM_H2);
		
		MolecularModel promModelProm_L_ACT = setupAcetyleneModel();
//		System.out.println(promModelProm_L_ACT);
		
		MolecularModel promModelProm_L_THF = setupTHFModel();
//		System.out.println(promModelProm_L_THF);	
		
		MolecularModel waterModelWaterSPCE = setupWaterModelSPCE(alw_SPCE, dHOw_SPCE);
		
		Integer keyModelMap = 0;
		molModelMap = new HashMap<Integer, MolecularModel>();
	
		molModelMap.put(1, waterModelWater);		
		molModelMap.put(nextKeyInMap(molModelMap), hydrogenModelProm_SM_H2);
		molModelMap.put(nextKeyInMap(molModelMap), promModelProm_L_ACT);
		molModelMap.put(nextKeyInMap(molModelMap), promModelProm_L_THF);
		molModelMap.put(nextKeyInMap(molModelMap), waterModelWaterSPCE);
		
//		ArrayList<Molecule> molecules = new ArrayList<Molecule>(numbMolsTotal);
		
		AtomVector molCenterOfMass;
		MolEulerAngles molEulerAngles;
		ArrayList<AtomVector> molLJSites;
		ArrayList<AtomVector> molCoulombSites;
		
		Molecule molecule;
		
		ArrayList<AtomVector> allAtoms = new ArrayList<AtomVector>(3*numbMolsTotal);		
		readXYZ_File("input_sH_clathrate_water_coords.xyz", allAtoms);
		
		int k = 0;
		
		for (int i=0; i<allAtoms.size(); i+=3){
			
			molCenterOfMass = new AtomVector(allAtoms.get(k+i));			
			molEulerAngles = new MolEulerAngles();
			
			AtomVector siteLJ;
			molLJSites = new ArrayList<AtomVector>(molModelMap.get(1).numbOfLJSites);
			
			siteLJ = new AtomVector(allAtoms.get(k+i));
			molLJSites.add(siteLJ);
			
			
			AtomVector siteCoulomb;
			molCoulombSites = new ArrayList<AtomVector>(molModelMap.get(1).numbOfCoulombSites);
			
			siteCoulomb = new AtomVector(allAtoms.get(k+i));
			molCoulombSites.add(siteCoulomb);
			
			siteCoulomb = new AtomVector(allAtoms.get(k+i+1));
			molCoulombSites.add(siteCoulomb);
			
			siteCoulomb = new AtomVector(allAtoms.get(k+i+2));
			molCoulombSites.add(siteCoulomb);
			
			molecule = new Molecule(1,
								    molModelMap.get(1).numbOfLJSites,
									molModelMap.get(1).numbOfCoulombSites,
									molCenterOfMass,
									molEulerAngles,
									molLJSites,
									molCoulombSites
								    );
			
			molecules.add(molecule);			
		}
		
		writeXYZ_VMD_File("output_initial_water_framework_VMD_1.xyz", molecules, molModelMap, 3*molecules.size());
		calculateWaterEulerXYZ(molecules, alw);
		writeXYZ_VMD_File("output_initial_water_framework_VMD_2.xyz", molecules, molModelMap, 3*molecules.size());
//		moleculesArrayCheck(molecules);
		
		allAtoms.clear();
		readXYZ_File("input_small_cages_centers_coords.xyz", allAtoms);
		readXYZ_File("input_medium_cages_centers_coords.xyz", allAtoms);
		
		
		ArrayList<Molecule> moleculesH2 = new ArrayList<Molecule>(numbHydrogenMols);		
		moleculesH2 = calculateGuestEulerXYZ(allAtoms,
				 							 molModelMap,
				 							 2, 
				 							 numbHydrogenMols
											);		
		int numbOfCoulombSites = molModelMap.get(2).numbOfCoulombSites;
		writeXYZ_VMD_File("output_initial_hyrogens_VMD.xyz", moleculesH2, molModelMap, numbOfCoulombSites*moleculesH2.size());
		
		
		allAtoms.clear();
		readXYZ_File("input_large_cages_centers_coords.xyz", allAtoms);		
		ArrayList<Molecule> moleculesACET = new ArrayList<Molecule>(numbLargeCages);		
		moleculesACET = calculateGuestEulerXYZ(allAtoms,
				 							   molModelMap,
				 							   3, 
				 							   numbLargeCages
											  );
		numbOfCoulombSites = molModelMap.get(3).numbOfCoulombSites;
		writeXYZ_VMD_File("output_initial_acetylene_VMD.xyz", moleculesACET, molModelMap, numbOfCoulombSites*moleculesACET.size());
		
		allAtoms.clear();
		readXYZ_File("input_large_cages_centers_coords.xyz", allAtoms);		
		ArrayList<Molecule> moleculesTHF = new ArrayList<Molecule>(numbLargeCages);		
		moleculesTHF = calculateGuestEulerXYZ(allAtoms,
				 							  molModelMap,
				 							  4, 
				 							  numbLargeCages
											  );
		numbOfCoulombSites = molModelMap.get(4).numbOfCoulombSites;
		writeXYZ_VMD_File("output_initial_THF_VMD.xyz", moleculesTHF, molModelMap, numbOfCoulombSites*moleculesTHF.size());

		// Total energy calculation
	
		volume = box.x * box.y * box.z;
		
		LinkedList<Double> listBox = new LinkedList<Double>();
		listBox.add(box.x);
		listBox.add(box.y);
		listBox.add(box.z);
		rcutSq = Collections.min(listBox);
		rcutSq = rcutSq*rcutSq;
		
		numbOfMolsActual = numbWaterMols;
		
		tableEnergyAllMoleciles =  new double [numbOfMolsActual][numbOfMolsActual];
		
		enTot = energyTotal(molModelMap, 
                		    molecules, 
                		    numbOfMolsActual, 
                		    box,
                		    rcutSq,
                		    tableEnergyAllMoleciles
						   );
		
		System.out.println("Initial Total Energy = " + enTot/(double)numbOfMolsActual);
		System.out.println("Initial box = " + box + "\n");
		
/*		for (int row=0; row < tableInteractionEnergy.length; row++){
			for (int col=0; col < tableInteractionEnergy[0].length; col++){
				System.out.println("i, j, table[][] = " + row + "\t" + col + "\t" + tableInteractionEnergy[row][col]);
			}
		}*/
		
	}
	
	protected final double volumeChange(final HashMap<Integer, MolecularModel> molModelMap, 
			    						ArrayList<Molecule> molecules, 
			    						final int numbOfMolsActual,
			    						AtomVector boxOld,
			    						final double rcutSq,
			    						double [][] tableEnergyAllMoleciles,
			    						double enTot
			   						   ){
		int start = -1;
	    int end = 1;
	    Random random = new Random();
	    
	    AtomVector boxNew = new AtomVector(boxOld.x + DL * RandomDouble(start, end, random),
	    								   boxOld.y + DL * RandomDouble(start, end, random),
	    								   boxOld.z + DL * RandomDouble(start, end, random)
	    								  );
	    
	    double volumeOld = boxOld.x * boxOld.y * boxOld.z;
	    double volumeNew = boxNew.x * boxNew.y * boxNew.z;
	    
	    int numbLJSitesMol;
	    int numbCoulSitesMol;   
	    
	    AtomVector siteLJ;
	    ArrayList<AtomVector> molLJSites;
	    
	    AtomVector siteCoulomb;
	    ArrayList<AtomVector> molCoulombSites;	    
	    
		ArrayList<AtomVector> moleculesCOMOld = new ArrayList<AtomVector>(numbOfMolsActual);		
		
		ArrayList<ArrayList<AtomVector>> moleculesLJSitesOld = new ArrayList<ArrayList<AtomVector>>(numbOfMolsActual);
		ArrayList<ArrayList<AtomVector>> moleculesCoulSitesOld = new ArrayList<ArrayList<AtomVector>>(numbOfMolsActual);
		
	    AtomVector molCenterOfMassOld;
	    AtomVector molCenterOfMassNew;	    
	    AtomVector diffCOM;
	    
	    for (int i=0; i < numbOfMolsActual; i++){

			molCenterOfMassOld = new AtomVector(molecules.get(i).molCenterOfMass);
			moleculesCOMOld.add(molCenterOfMassOld);
			
			molCenterOfMassNew = new AtomVector (molCenterOfMassOld.x * boxNew.x / boxOld.x,
												 molCenterOfMassOld.y * boxNew.y / boxOld.y,
												 molCenterOfMassOld.z * boxNew.z / boxOld.z
												);
			molecules.get(i).molCenterOfMass = molCenterOfMassNew.AtomVectorCopy();
			
			diffCOM = new AtomVector(molCenterOfMassNew.AtomVectorSubtract(molCenterOfMassOld));
			
	        numbLJSitesMol = molecules.get(i).numbOfLJSites;
	        numbCoulSitesMol = molecules.get(i).numbOfCoulombSites; 
			
			molLJSites = new ArrayList<AtomVector>();
	        for (int j=0; j < numbLJSitesMol; j++){
	        	
	        	siteLJ = new AtomVector(molecules.get(i).molLJSites.get(j));
	        	molLJSites.add(siteLJ);	    
	        }
	        moleculesLJSitesOld.add(molLJSites);
	        
	        for (int j=0; j < numbLJSitesMol; j++){
	        	
	        	siteLJ = new AtomVector(molecules.get(i).molLJSites.get(j));
	        	molecules.get(i).molLJSites.set(j, siteLJ.AtomVectorAdd(diffCOM));	        	
	        }	        
	        
	        molCoulombSites  = new ArrayList<AtomVector>();
	        for (int j=0; j < numbCoulSitesMol; j++){
	        	
	        	siteCoulomb = new AtomVector(molecules.get(i).molCoulombSites.get(j));
	        	molCoulombSites.add(siteCoulomb);
	        }
	        moleculesCoulSitesOld.add(molCoulombSites);       
	        
	        for (int j=0; j < numbCoulSitesMol; j++){
	        	
	        	siteCoulomb = new AtomVector(molecules.get(i).molCoulombSites.get(j));
	        	molecules.get(i).molCoulombSites.set(j, siteCoulomb.AtomVectorAdd(diffCOM));
	        }	        
	    }
	        
	    double Ua = enTot;
	    
	    double [][] tableEnergyAllMolecilesNew = new double [numbOfMolsActual][numbOfMolsActual];
	        
		double Ub = energyTotal(molModelMap, 
        	    				molecules, 
        	    				numbOfMolsActual, 
        	    				boxNew,
        	    				rcutSq,
        	    				tableEnergyAllMolecilesNew
							   );
			
		double boltz1 = numbOfMolsActual * Math.log(volumeNew/volumeOld);
		double boltz2 = -0.0143836 * pressure * (volumeNew - volumeOld) * invTemperature * KBR;
		double boltz3 = -invTemperature * (Ub - Ua) * KBR;

	    double prob = Math.exp(boltz1 + boltz2 + boltz3);
	        
		start = 0;
		end = 1;		
		double ran = RandomDouble(start, end, random);   
		    
		if (ran < prob){
		    	
		   	enTot = Ub;
		    	
		   	for(int i=0; i < numbOfMolsActual; i++){		    		
		   		for(int j=0; j < numbOfMolsActual; j++){
		   			
		   			tableEnergyAllMoleciles[i][j] = tableEnergyAllMolecilesNew[i][j];		   			
		   		}		    		
		   	}
		    	
//	    	boxOld = boxNew;  	
		   	boxOld.x = boxNew.x;
		   	boxOld.y = boxNew.y;
		   	boxOld.z = boxNew.z;
		   	
	    	volume = volumeNew;
		    	
		    } else {
		    	
		    	for (int i=0; i < numbOfMolsActual; i++){
		    		
			        numbLJSitesMol = molecules.get(i).numbOfLJSites;
			        numbCoulSitesMol = molecules.get(i).numbOfCoulombSites; 
			        					
			        for (int j=0; j < numbLJSitesMol; j++){
			        	molecules.get(i).molLJSites.set(j, moleculesLJSitesOld.get(i).get(j));
			        }	        
			        
			        for (int j = 0; j < numbCoulSitesMol; j++) {
			        	siteCoulomb = new AtomVector(molecules.get(i).molCoulombSites.get(j));

			        }
		    	}  	
		    }	
		return enTot;
	}
	
	
	protected final double moleculeRotation(final int imol,
				 						    final HashMap<Integer, MolecularModel> molModelMap, 
				 						    ArrayList<Molecule> molecules, 
				 						    final int numbOfMolsActual,
				 						    final AtomVector box,
				 						    final double rcutSq,
				 						    double [][] tableEnergyAllMoleciles,
				 						    double enTot
				 						   ){
		
		Integer keyMapMolI;
		
        keyMapMolI = molecules.get(imol).molModelKey;
        int numbLJSitesMolI = molModelMap.get(keyMapMolI).numbOfLJSites;
        int numbCoulSitesMolI = molModelMap.get(keyMapMolI).numbOfCoulombSites; 
		
		ArrayList<AtomVector> molLJSitesOld = new ArrayList<AtomVector>(numbLJSitesMolI);	
		molLJSitesOld.addAll(molecules.get(imol).molLJSites);		
/*		for (int iatom = 0; iatom < numbLJSitesMolI; iatom++) {			
			molLJSitesOld.add(molecules.get(imol).molLJSites.get(iatom));
		}*/
		
		ArrayList<AtomVector> molCoulombSitesOld = new ArrayList<AtomVector>(numbCoulSitesMolI);
		molCoulombSitesOld.addAll(molecules.get(imol).molCoulombSites);
/*		for (int iatom = 0; iatom < numbCoulSitesMolI; iatom++) {
			molCoulombSitesOld.add(molecules.get(imol).molCoulombSites.get(iatom));
		}*/
		
	    double alphaTemp = molecules.get(imol).molEulerAngles.alpha;
	    double betaTemp  = molecules.get(imol).molEulerAngles.beta;
	    double gammaTemp = molecules.get(imol).molEulerAngles.gamma;
	    
	    double alphaOld = alphaTemp;
	    double betaOld  = betaTemp;
	    double gammaOld = gammaTemp;
	    
		int start = -1;
	    int end = 1;
	    Random random = new Random();

	    alphaTemp += MAXANG * RandomDouble(start, end, random);
	    betaTemp  += MAXCOS * RandomDouble(start, end, random);
	    gammaTemp += MAXANG * RandomDouble(start, end, random);
	    
		// keep alphaTemp and gammaTemp angles in [-PI, PI]
		// keep betaTemp = dcos(betaTemp) in [-1,1]
	    alphaTemp = alphaTemp - Math.rint(alphaTemp / TWOPI) * TWOPI;
	    betaTemp  = betaTemp  - Math.rint(betaTemp / 2.0) * 2.0;
	    gammaTemp = gammaTemp - Math.rint(gammaTemp / TWOPI) * TWOPI;
	    
		double cosAlpha, cosBeta, cosGamma;
		double sinAlpha, sinBeta, sinGamma;
	    
        cosAlpha = Math.cos(alphaTemp);
        cosBeta  = betaTemp;
        cosGamma = Math.cos(gammaTemp);

        sinAlpha = Math.sin(alphaTemp);
        sinBeta  = Math.sqrt(1.0-cosBeta*cosBeta);
        sinGamma = Math.sin(gammaTemp);
        
		double A11, A12, A13;
		double A21, A22, A23;
		double A31, A32, A33;
        
        A11 =  cosAlpha*cosGamma - sinAlpha*cosBeta*sinGamma;
        A12 = -cosAlpha*sinGamma - sinAlpha*cosBeta*cosGamma;
        A13 =  sinBeta*sinAlpha;

        A21 =  sinAlpha*cosGamma + cosAlpha*cosBeta*sinGamma;
        A22 = -sinAlpha*sinGamma + cosAlpha*cosBeta*cosGamma;
        A23 = -sinBeta*cosAlpha;

        A31 =  sinBeta*sinGamma;
        A32 =  sinBeta*cosGamma;
        A33 =  cosBeta;
        
        AtomVector siteLJ, siteCoulomb, molCOM;
        AtomVector fixedAtom;
        
        molCOM = molecules.get(imol).molCenterOfMass;
        
        for (int j=0; j < numbLJSitesMolI; j++){
        	
        	fixedAtom = molModelMap.get(keyMapMolI).sitesLJ.get(j).siteLJ;
        	
        	siteLJ = new AtomVector(); 
        	
			siteLJ.x = A11*fixedAtom.x + A12*fixedAtom.y + A13*fixedAtom.z + molCOM.x;  
			siteLJ.y = A21*fixedAtom.x + A22*fixedAtom.y + A23*fixedAtom.z + molCOM.y;
			siteLJ.z = A31*fixedAtom.x + A32*fixedAtom.y + A33*fixedAtom.z + molCOM.z;
			
			molecules.get(imol).molLJSites.set(j, siteLJ);        	
        }
                       
        for (int j=0; j < numbCoulSitesMolI; j++){
        	
        	fixedAtom = molModelMap.get(keyMapMolI).sitesCoulomb.get(j).siteCoulomb;
        	
        	siteCoulomb = new AtomVector(); 
        	
        	siteCoulomb.x = A11*fixedAtom.x + A12*fixedAtom.y + A13*fixedAtom.z + molCOM.x;  
        	siteCoulomb.y = A21*fixedAtom.x + A22*fixedAtom.y + A23*fixedAtom.z + molCOM.y;
        	siteCoulomb.z = A31*fixedAtom.x + A32*fixedAtom.y + A33*fixedAtom.z + molCOM.z;
			
			molecules.get(imol).molCoulombSites.set(j, siteCoulomb);        	
        }	
        
		double  Ua = 0.0, Ub = 0.0;
		
		for (int jmol=0; jmol < numbOfMolsActual; jmol++){
			
			if (jmol > imol){
				Ua += tableEnergyAllMoleciles[imol][jmol];
			}else if(jmol < imol){
				Ua += tableEnergyAllMoleciles[jmol][imol];
			}			
		}
		
		double [] tableEnergyOneMolecile = new double[numbOfMolsActual];
		
        Ub = energyMolecule(imol,
				  			molModelMap, 
				  			molecules, 
				  			numbOfMolsActual, 
				  			box,
				  			rcutSq,
				  			tableEnergyOneMolecile
				 			);
        
		start = 0;
	    end = 1;		
	    double ran = RandomDouble(start, end, random);
	    
	    double prob = Math.exp(invTemperature * (Ua - Ub) * TOCAL);
	    
	    if (ran < prob){
	    	
	        molecules.get(imol).molEulerAngles.alpha = alphaOld; 
		    molecules.get(imol).molEulerAngles.beta  = betaOld;  
		    molecules.get(imol).molEulerAngles.gamma = gammaOld;   
	    	
			for (int jmol=0; jmol < numbOfMolsActual; jmol++){
				
				if (jmol > imol){
					tableEnergyAllMoleciles[imol][jmol] = tableEnergyOneMolecile[jmol];
				}else if(jmol < imol){
					tableEnergyAllMoleciles[jmol][imol] = tableEnergyOneMolecile[jmol];
				}			
			}			
			enTot += -Ua + Ub;			
	    } else {
	    	
			for (int iatom = 0; iatom < numbLJSitesMolI; iatom++) {			
				molecules.get(imol).molLJSites.set(iatom, molLJSitesOld.get(iatom));
			}
			
 			for (int iatom = 0; iatom < numbCoulSitesMolI; iatom++) {		
				molecules.get(imol).molCoulombSites.set(iatom, molCoulombSitesOld.get(iatom));
			}
	    }
        
		return enTot;
	}	
	
	protected final double moleculeTranslation(final int imol,
			  								   final HashMap<Integer, MolecularModel> molModelMap, 
			  								   ArrayList<Molecule> molecules, 
			  								   final int numbOfMolsActual,
			  								   final AtomVector box,
			  								   final double rcutSq,
			  								   double [][] tableEnergyAllMoleciles,
			  								   double enTot
											  ){
		
		Integer keyMapMolI, keyMapMolJ;
		int numbLJSitesMolI, numbLJSitesMolJ;
		int numbCoulSitesMolI, numbCoulSitesMolJ;
		
		AtomVector molCenterOfMassOld = new AtomVector();	
		molCenterOfMassOld = molecules.get(imol).molCenterOfMass;
		
		numbLJSitesMolI = molecules.get(imol).numbOfLJSites;
		numbCoulSitesMolI = molecules.get(imol).numbOfCoulombSites;
		
		ArrayList<AtomVector> molLJSitesOld = new ArrayList<AtomVector>(numbLJSitesMolI);	
		molLJSitesOld.addAll(molecules.get(imol).molLJSites);		
/*		for (int iatom = 0; iatom < numbLJSitesMolI; iatom++) {			
			molLJSitesOld.add(molecules.get(imol).molLJSites.get(iatom));
		}*/
		
		ArrayList<AtomVector> molCoulombSitesOld = new ArrayList<AtomVector>(numbCoulSitesMolI);
		molCoulombSitesOld.addAll(molecules.get(imol).molCoulombSites);
/*		for (int iatom = 0; iatom < numbCoulSitesMolI; iatom++) {
			molCoulombSitesOld.add(molecules.get(imol).molCoulombSites.get(iatom));
		}*/

		int start = -1;
	    int end = 1;
	    Random random = new Random();
	    
	    AtomVector vecRandomTransl = new AtomVector(MAXTRANS * RandomDouble(start, end, random),
	    										    MAXTRANS * RandomDouble(start, end, random),
	    										    MAXTRANS * RandomDouble(start, end, random)
	    										   );
	    
	    molecules.get(imol).molCenterOfMass.AtomVectorAdd(vecRandomTransl);
	    
		for (int iatom = 0; iatom < numbLJSitesMolI; iatom++) {			
			molecules.get(imol).molLJSites.set(iatom, 
					                           molecules.get(imol).molLJSites.get(iatom).AtomVectorAdd(vecRandomTransl)
					                          );
		}
		
		for (int iatom = 0; iatom < numbCoulSitesMolI; iatom++) {		
			molecules.get(imol).molCoulombSites.set(iatom,
					                                molecules.get(imol).molCoulombSites.get(iatom).AtomVectorAdd(vecRandomTransl)
					                               );
		}
		
		double  Ua = 0.0, Ub = 0.0;
		
		for (int jmol=0; jmol < numbOfMolsActual; jmol++){
			
			if (jmol > imol){
				Ua += tableEnergyAllMoleciles[imol][jmol];
			}else if(jmol < imol){
				Ua += tableEnergyAllMoleciles[jmol][imol];
			}			
		}
		
		double [] tableEnergyOneMolecile = new double[numbOfMolsActual];
		
        Ub = energyMolecule(imol,
				  			molModelMap, 
				  			molecules, 
				  			numbOfMolsActual, 
				  			box,
				  			rcutSq,
				  			tableEnergyOneMolecile
				 			);
        
		start = 0;
	    end = 1;		
	    double ran = RandomDouble(start, end, random);
	    
	    double prob = Math.exp(invTemperature * (Ua - Ub) * TOCAL);
	    
	    if (ran < prob){
	    	
			for (int jmol=0; jmol < numbOfMolsActual; jmol++){
				
				if (jmol > imol){
					tableEnergyAllMoleciles[imol][jmol] = tableEnergyOneMolecile[jmol];
				}else if(jmol < imol){
					tableEnergyAllMoleciles[jmol][imol] = tableEnergyOneMolecile[jmol];
				}			
			}			
			enTot += -Ua + Ub;			
	    } else {
	    	
			molecules.get(imol).molCenterOfMass = molCenterOfMassOld;
			
			for (int iatom = 0; iatom < numbLJSitesMolI; iatom++) {			
				molecules.get(imol).molLJSites.set(iatom, molLJSitesOld.get(iatom));
			}
		
 			for (int iatom = 0; iatom < numbCoulSitesMolI; iatom++) {		
				molecules.get(imol).molCoulombSites.set(iatom, molCoulombSitesOld.get(iatom));
			}
	    }
	    return enTot;
	}
	
	protected final double energyMolecule(final int imol,
										  final HashMap<Integer, MolecularModel> molModelMap, 
										  final ArrayList<Molecule> molecules, 
										  final int numbOfMolsActual, 
										  final AtomVector box,
										  final double rcutSq,
										  double [] tableEnergyOneMolecile
										 ){

		double enMolU = 0.0;
		double enLJ, enCoul;
		Integer keyMapMolI, keyMapMolJ;
		int numbLJSitesMolI, numbLJSitesMolJ;
		int numbCoulSitesMolI, numbCoulSitesMolJ;
		double r, rSq, r6, r12;
		double epsilon, epsilonI, epsilonJ, sigma, sigmaI, sigmaJ, sigma6, sigma12, attr, repul;
		double chargeI, chargeJ;

		AtomVector atomI = new AtomVector();
		AtomVector atomJ = new AtomVector();
		AtomVector vecIJ = new AtomVector();
		AtomVector pbc   = new AtomVector();

		enMolU = 0.0;

		keyMapMolI = molecules.get(imol).molModelKey;
		numbLJSitesMolI = molModelMap.get(keyMapMolI).numbOfLJSites;
		numbCoulSitesMolI = molModelMap.get(keyMapMolI).numbOfCoulombSites;
		
		for (int jmol = 0; jmol < numbOfMolsActual; jmol++) {
			
			if ( jmol == imol ){
				
				tableEnergyOneMolecile[jmol] = 0.0;
				
			}else{

				tableEnergyOneMolecile[jmol] = 0.0;

				keyMapMolJ = molecules.get(jmol).molModelKey;
				numbLJSitesMolJ = molecules.get(jmol).numbOfLJSites;

				for (int iatom = 0; iatom < numbLJSitesMolI; iatom++) {
					for (int jatom = 0; jatom < numbLJSitesMolJ; jatom++) {

						atomI = molecules.get(imol).molLJSites.get(iatom);
						atomJ = molecules.get(jmol).molLJSites.get(jatom);

						vecIJ = atomI.AtomVectorSubtract(atomJ);

						pbc.x = box.x * Math.rint(vecIJ.x / box.x);
						pbc.y = box.y * Math.rint(vecIJ.y / box.y);
						pbc.z = box.z * Math.rint(vecIJ.z / box.z);

						vecIJ = vecIJ.AtomVectorSubtract(pbc);
						rSq = vecIJ.AtomVectorLengthSquared();

						if (rSq < rcutSq) {

							r6 = rSq * rSq * rSq;
							r12 = r6 * r6;

							sigmaI = molModelMap.get(keyMapMolI).sitesLJ.get(iatom).sigma;
							sigmaJ = molModelMap.get(keyMapMolJ).sitesLJ.get(jatom).sigma;

							sigma = 0.5 * (sigmaI + sigmaJ);

							sigma6 = Math.pow(sigma, 6.0);
							sigma12 = sigma6 * sigma6;

							epsilonI = molModelMap.get(keyMapMolI).sitesLJ.get(iatom).epsilon;
							epsilonJ = molModelMap.get(keyMapMolJ).sitesLJ.get(jatom).epsilon;

							epsilon = Math.sqrt(epsilonI * epsilonJ);

							attr = sigma6 / r6;
							repul = sigma12 / r12;

							enLJ = 4.0 * epsilon * (repul - attr);

							tableEnergyOneMolecile[jmol] += enLJ;
							enMolU += enLJ;
						}
					}
				}

				numbCoulSitesMolJ = molModelMap.get(keyMapMolJ).numbOfCoulombSites;

				for (int iatom = 0; iatom < numbCoulSitesMolI; iatom++) {
					for (int jatom = 0; jatom < numbCoulSitesMolJ; jatom++) {

						atomI = molecules.get(imol).molCoulombSites.get(iatom);
						atomJ = molecules.get(jmol).molCoulombSites.get(jatom);

						vecIJ = atomI.AtomVectorSubtract(atomJ);

						pbc.x = box.x * Math.rint(vecIJ.x / box.x);
						pbc.y = box.y * Math.rint(vecIJ.y / box.y);
						pbc.z = box.z * Math.rint(vecIJ.z / box.z);

						vecIJ = vecIJ.AtomVectorSubtract(pbc);
						rSq = vecIJ.AtomVectorLengthSquared();

						if (rSq < rcutSq) {

							r = Math.sqrt(rSq);

							chargeI = molModelMap.get(keyMapMolI).sitesCoulomb.get(iatom).charge;
							chargeJ = molModelMap.get(keyMapMolJ).sitesCoulomb.get(jatom).charge;

							enCoul = TOCAL * chargeI * chargeJ / r;

							tableEnergyOneMolecile[jmol] += enCoul;
							enMolU += enCoul;
						}
					}
				}
			}	
		}
		return enMolU;
	}
	
	protected final double energyTotal(final HashMap<Integer, MolecularModel> molModelMap, 
			                           final ArrayList<Molecule> molecules, 
			                           final int numbOfMolsActual, 
			                           final AtomVector box,
			                           final double rcutSq,
			                           double [][] tableEnergyAllMoleciles
			                          ){
		
		double enTot = 0.0;
		double enLJ, enCoul;
		Integer keyMapMolI, keyMapMolJ;
		int numbLJSitesMolI, numbLJSitesMolJ;
		int numbCoulSitesMolI, numbCoulSitesMolJ;
		double r, rSq, r6, r12;
		double epsilon, epsilonI, epsilonJ, sigma, sigmaI, sigmaJ, sigma6, sigma12, attr, repul;
		double chargeI, chargeJ;
		
		AtomVector atomI = new AtomVector();
		AtomVector atomJ = new AtomVector();
		AtomVector vecIJ = new AtomVector();
		AtomVector pbc   = new AtomVector();
		
		enTot = 0.0;
		
		for (int imol=0; imol<numbOfMolsActual - 1; imol++){
			for (int jmol=imol+1; jmol<numbOfMolsActual; jmol++){
				
				tableEnergyAllMoleciles[imol][jmol] = 0.0;
				
				keyMapMolI = molecules.get(imol).molModelKey;
				numbLJSitesMolI = molModelMap.get(keyMapMolI).numbOfLJSites;
				
				keyMapMolJ = molecules.get(jmol).molModelKey;
				numbLJSitesMolJ = molModelMap.get(keyMapMolJ).numbOfLJSites;
				
				for (int iatom=0; iatom < numbLJSitesMolI; iatom++){
					for (int jatom=0; jatom < numbLJSitesMolJ; jatom++){
						
						atomI =  molecules.get(imol).molLJSites.get(iatom);    
						atomJ =  molecules.get(jmol).molLJSites.get(jatom);
				
						vecIJ = atomI.AtomVectorSubtract(atomJ);
						
						pbc.x = box.x*Math.rint(vecIJ.x/box.x);
						pbc.y = box.y*Math.rint(vecIJ.y/box.y);
						pbc.z = box.z*Math.rint(vecIJ.z/box.z);
						
						vecIJ = vecIJ.AtomVectorSubtract(pbc);						
						rSq = vecIJ.AtomVectorLengthSquared();
						
						if (rSq < rcutSq){
							
							r6 = rSq*rSq*rSq;
							r12 = r6*r6;
							
							sigmaI = molModelMap.get(keyMapMolI).sitesLJ.get(iatom).sigma;
							sigmaJ = molModelMap.get(keyMapMolJ).sitesLJ.get(jatom).sigma;
							
							sigma = 0.5*(sigmaI + sigmaJ);
							
		                    sigma6 = Math.pow(sigma, 6.0);
		                    sigma12 = sigma6*sigma6;
		                    
		                    epsilonI = molModelMap.get(keyMapMolI).sitesLJ.get(iatom).epsilon; 
		                    epsilonJ = molModelMap.get(keyMapMolJ).sitesLJ.get(jatom).epsilon;		
		                    
		                    epsilon = Math.sqrt(epsilonI * epsilonJ);
		                    
		                    attr = sigma6/r6;
		                    repul = sigma12/r12;
		                    
		                    enLJ = 4.0 * epsilon * (repul - attr);
							
		                    tableEnergyAllMoleciles[imol][jmol] += enLJ;		                    
		                    enTot += enLJ;
						}	
					}
				}
				
				numbCoulSitesMolI = molModelMap.get(keyMapMolI).numbOfCoulombSites;
				numbCoulSitesMolJ = molModelMap.get(keyMapMolJ).numbOfCoulombSites;
				
				for (int iatom=0; iatom < numbCoulSitesMolI; iatom++){
					for (int jatom=0; jatom < numbCoulSitesMolJ; jatom++){
						
						atomI =  molecules.get(imol).molCoulombSites.get(iatom);    
						atomJ =  molecules.get(jmol).molCoulombSites.get(jatom);
						
						vecIJ = atomI.AtomVectorSubtract(atomJ);
						
						pbc.x = box.x*Math.rint(vecIJ.x/box.x);
						pbc.y = box.y*Math.rint(vecIJ.y/box.y);
						pbc.z = box.z*Math.rint(vecIJ.z/box.z);
						
						vecIJ = vecIJ.AtomVectorSubtract(pbc);						
						rSq = vecIJ.AtomVectorLengthSquared();
						
						if (rSq < rcutSq){
							
							r = Math.sqrt(rSq);
							
							chargeI = molModelMap.get(keyMapMolI).sitesCoulomb.get(iatom).charge;
							chargeJ = molModelMap.get(keyMapMolJ).sitesCoulomb.get(jatom).charge;
							
							enCoul = TOCAL * chargeI * chargeJ / r;
							
							tableEnergyAllMoleciles[imol][jmol] += enCoul;		                    
		                    enTot += enCoul;							
						}						
					}
				}				
			}			
		}		
		return enTot;
	}
	
	
	private static int RandomInteger(int aStart, int aEnd, Random aRandom) {
		if (aStart > aEnd) {
			throw new IllegalArgumentException("Start cannot exceed End.");
		}
		// get the range, casting to long to avoid overflow problems
		long range = (long) aEnd - (long) aStart + 1;
		// compute a fraction of the range, 0 <= frac < range
		long fraction = (long) (range * aRandom.nextDouble());
		int randomNumber = (int) (fraction + aStart);
		return randomNumber;
	}
	
	private static double RandomDouble(double aStart, double aEnd, Random aRandom) {
		if (aStart > aEnd) {
			throw new IllegalArgumentException("Start cannot exceed End.");
		}
		double range = aEnd - aStart;
		// compute a fraction of the range, 0 <= frac < range
		double fraction = (range * aRandom.nextDouble());
		double randomNumber = fraction + aStart;
		return randomNumber;
	}
 

	protected ArrayList<Molecule> calculateGuestEulerXYZ(final ArrayList<AtomVector> cages,
														 final HashMap<Integer, MolecularModel> molModelMap,
														 final Integer keyMolModel, 
														 final int numbOfGuests
														){
		
		double eAlpha, eBeta, eGamma;
		
		double cosAlpha, cosBeta, cosGamma;
		double sinAlpha, sinBeta, sinGamma;
		
		double A11, A12, A13;
		double A21, A22, A23;
		double A31, A32, A33;
		
		AtomVector molCenterOfMass;
		MolEulerAngles molEulerAngles;
		ArrayList<AtomVector> molLJSites;
		ArrayList<AtomVector> molCoulombSites;
		Molecule molecule;
		
		ArrayList<Molecule> mols = new ArrayList<Molecule>(numbOfGuests);
		
		int start = -1;
	    int end = 1;
	    Random random = new Random();
	    
	    for (int i=0; i<numbOfGuests; i++){
	    	
	    	molCenterOfMass = new AtomVector(cages.get(i));   // center of mass of molecule is put in the center of a cage
		
	    	// generate random orientation of molecule
			eAlpha = TWOPI * RandomDouble(start, end, random);              // Euler angle alpha
			eBeta = Math.cos((Math.PI) * RandomDouble(start, end, random)); // Euler angle beta
			eGamma = TWOPI * RandomDouble(start, end, random);              // Euler angle gamma

			// keep eAlpha and eGamma angles in [-PI, PI]
			// keep eBeta(i) = dcos(e_beta) in [-1,1]
			eAlpha = eAlpha - Math.rint(eAlpha / TWOPI) * TWOPI;
			eBeta  = eBeta - Math.rint(eBeta / 2.0) * 2.0;
			eGamma = eGamma - Math.rint(eGamma / TWOPI) * TWOPI;
			
			molEulerAngles = new MolEulerAngles();

			molEulerAngles.alpha = eAlpha;
			molEulerAngles.beta  = eBeta;
			molEulerAngles.gamma = eGamma;
			
            cosAlpha = Math.cos(eAlpha);
            cosBeta  = eBeta;
            cosGamma = Math.cos(eGamma);

            sinAlpha = Math.sin(eAlpha);
            sinBeta  = Math.sqrt(1.0-cosBeta*cosBeta);
            sinGamma = Math.sin(eGamma);
			
            A11 =  cosAlpha*cosGamma - sinAlpha*cosBeta*sinGamma;
            A12 = -cosAlpha*sinGamma - sinAlpha*cosBeta*cosGamma;
            A13 =  sinBeta*sinAlpha;

            A21 =  sinAlpha*cosGamma + cosAlpha*cosBeta*sinGamma;
            A22 = -sinAlpha*sinGamma + cosAlpha*cosBeta*cosGamma;
            A23 = -sinBeta*cosAlpha;

            A31 =  sinBeta*sinGamma;
            A32 =  sinBeta*cosGamma;
            A33 =  cosBeta;

            int numbOfLJSites = molModelMap.get(keyMolModel).numbOfLJSites;
            molLJSites = new ArrayList<AtomVector>(numbOfLJSites);
            AtomVector siteLJ;
            AtomVector fixedAtom;
            
            for (int j=0; j<numbOfLJSites; j++){
            	
            	fixedAtom = molModelMap.get(keyMolModel).sitesLJ.get(j).siteLJ;
            	
            	siteLJ = new AtomVector(); 
            	
    			siteLJ.x = A11*fixedAtom.x + A12*fixedAtom.y + A13*fixedAtom.z + cages.get(i).x;  
    			siteLJ.y = A21*fixedAtom.x + A22*fixedAtom.y + A23*fixedAtom.z + cages.get(i).y;
    			siteLJ.z = A31*fixedAtom.x + A32*fixedAtom.y + A33*fixedAtom.z + cages.get(i).z;
            	
    			molLJSites.add(siteLJ);
            	
            }
            
           
            int numbOfCoulombSites = molModelMap.get(keyMolModel).numbOfCoulombSites;  
            molCoulombSites = new ArrayList<AtomVector>(numbOfCoulombSites);
            AtomVector siteCoulomb;
            
            for (int j=0; j<numbOfCoulombSites; j++){
            	
            	fixedAtom = molModelMap.get(keyMolModel).sitesCoulomb.get(j).siteCoulomb;
            	
            	siteCoulomb = new AtomVector();
            	
    			siteCoulomb.x = A11*fixedAtom.x + A12*fixedAtom.y + A13*fixedAtom.z + cages.get(i).x;  
    			siteCoulomb.y = A21*fixedAtom.x + A22*fixedAtom.y + A23*fixedAtom.z + cages.get(i).y;
    			siteCoulomb.z = A31*fixedAtom.x + A32*fixedAtom.y + A33*fixedAtom.z + cages.get(i).z;
            	
    			molCoulombSites.add(siteCoulomb);
            	
            }
            
			molecule = new Molecule(keyMolModel,
				    				numbOfLJSites,
				    				numbOfCoulombSites,
				    				molCenterOfMass,
				    				molEulerAngles,
				    				molLJSites,
				    				molCoulombSites
				    			   );
			mols.add(molecule);		
	    }			
		
		return mols;
	}
	
	protected void calculateWaterEulerXYZ(final ArrayList<Molecule> molecules,
			                              final double alw
			                             ){
		
		AtomVector u,v, X, Y, Z;
		double norm;
		AtomVector vec1, vec2;
		AtomVector vecOxygen, vecHydrogen1, vecHydrogen2;
		
		double eAlpha, eBeta, eGamma;
		
		double cosAlpha, cosBeta, cosGamma;
		double sinAlpha, sinBeta, sinGamma;
		
		double A11, A12, A13;
		double A21, A22, A23;
		double A31, A32, A33;
		
//		System.exit(0);
		
		for (int i=0; i<numbWaterMols; i++){
			
			vecOxygen    = molecules.get(i).molCoulombSites.get(0);
			vecHydrogen1 = molecules.get(i).molCoulombSites.get(1);
			vecHydrogen2 = molecules.get(i).molCoulombSites.get(2);
		
			u = vecHydrogen1.AtomVectorSubtract(vecOxygen);
			v = vecHydrogen2.AtomVectorSubtract(vecOxygen);	
			
			Z = u.crossProduct(v);
			norm = Z.AtomVectorLength();
			Z = Z.AtomVectorDevideByScalar(norm);
			
			vec1 = vecOxygen;
			vec2 = Z.AtomVectorAdd(vecOxygen); 
			
			double theta = -(90.0 - alw/2.0)*Math.PI/180.0;  // angle of rotation about axis Z, in radians			
			X = rotatePoint(vecHydrogen1, theta, vec1, vec2);
			
			theta = -alw/2.0*Math.PI/180.0;  // angle of rotation about axis Z, in radians
			Y = rotatePoint(vecHydrogen2, theta, vec1, vec2);
			
			X.z = X.z - vecOxygen.z;
			Y.z = Y.z - vecOxygen.z;
			
//			get the Euler Angles of a given frame
//			http://en.wikipedia.org/wiki/Euler_angles#Geometric_derivation			
			eAlpha = Math.atan2(Z.x,-Z.y);                		  // Euler angle alpha
			eBeta  = Math.atan2(Math.sqrt(Z.x*Z.x+Z.y*Z.y),Z.z);  // Euler angle beta
			eGamma = Math.atan2(X.z,Y.z);                 		  // Euler angle gamma

			eBeta = Math.cos(eBeta);
			
//			keep eAlpha and eGamma angles in [-PI, PI]
//			keep eBeta(i) = dcos(e_beta) in [-1,1]			
			eAlpha = eAlpha - Math.rint(eAlpha/TWOPI)*TWOPI;
			eBeta  = eBeta  - Math.rint(eBeta/2.0)*2.0;
			eGamma = eGamma - Math.rint(eGamma/TWOPI)*TWOPI;
					            
			molecules.get(i).molEulerAngles.alpha = eAlpha;
			molecules.get(i).molEulerAngles.beta  = eBeta;
			molecules.get(i).molEulerAngles.gamma = eGamma;
			
            cosAlpha = Math.cos(eAlpha);
            cosBeta  = eBeta;
            cosGamma = Math.cos(eGamma);

            sinAlpha = Math.sin(eAlpha);
            sinBeta  = Math.sqrt(1.0-cosBeta*cosBeta);
            sinGamma = Math.sin(eGamma);
			
            A11 =  cosAlpha*cosGamma - sinAlpha*cosBeta*sinGamma;
            A12 = -cosAlpha*sinGamma - sinAlpha*cosBeta*cosGamma;
            A13 =  sinBeta*sinAlpha;

            A21 =  sinAlpha*cosGamma + cosAlpha*cosBeta*sinGamma;
            A22 = -sinAlpha*sinGamma + cosAlpha*cosBeta*cosGamma;
            A23 = -sinBeta*cosAlpha;

            A31 =  sinBeta*sinGamma;
            A32 =  sinBeta*cosGamma;
            A33 =  cosBeta;
            
            Integer keyMap = molecules.get(i).molModelKey;
            AtomVector fixedHydrogen1 =  molModelMap.get(keyMap).sitesCoulomb.get(1).siteCoulomb;
            AtomVector fixedHydrogen2 =  molModelMap.get(keyMap).sitesCoulomb.get(2).siteCoulomb;
            
            // Hydrohen1
			molecules.get(i).molCoulombSites.get(1).x = A11*fixedHydrogen1.x + A12*fixedHydrogen1.y + A13*fixedHydrogen1.z + vecOxygen.x;  
			molecules.get(i).molCoulombSites.get(1).y = A21*fixedHydrogen1.x + A22*fixedHydrogen1.y + A23*fixedHydrogen1.z + vecOxygen.y;
			molecules.get(i).molCoulombSites.get(1).z = A31*fixedHydrogen1.x + A32*fixedHydrogen1.y + A33*fixedHydrogen1.z + vecOxygen.z;
			
			// Hydrohen2
			molecules.get(i).molCoulombSites.get(2).x = A11*fixedHydrogen2.x + A12*fixedHydrogen2.y + A13*fixedHydrogen2.z + vecOxygen.x;
			molecules.get(i).molCoulombSites.get(2).y = A21*fixedHydrogen2.x + A22*fixedHydrogen2.y + A23*fixedHydrogen2.z + vecOxygen.y;
			molecules.get(i).molCoulombSites.get(2).z = A31*fixedHydrogen2.x + A32*fixedHydrogen2.y + A33*fixedHydrogen2.z + vecOxygen.z;
		
		}		
	}

	
	protected void moleculesArrayCheck(final ArrayList<Molecule> mols){
		
		for (int i = 0; i < mols.size(); i++) {
			System.out.println(i + "\t" + mols.get(i).molEulerAngles);
			
		}		
	}
	
/*	     Rotation of a point vec in 3 dimensional space by angle
	     theta about an arbitrary axes defined by a line between two points vec1 and vec2 */
	
	protected AtomVector rotatePoint(final AtomVector vec0, final double theta, final AtomVector vec1, final AtomVector vec2){
		
		AtomVector vec, u, vectorResult;
        double  norm,c,s;
        double qx, qy, qz;
        double M11, M12, M13;
        double M21, M22, M23;
        double M31, M32, M33;

        u = vec2.AtomVectorSubtract(vec1); 
        
        norm = u.AtomVectorLength();        		
        u = u.AtomVectorDevideByScalar(norm); 
        
        c = Math.cos(theta);
        s = Math.sin(theta);
        
        vec = vec0;
        vec = vec.AtomVectorSubtract(vec1);
 
        M11 = u.x*u.x*(1-c)+c;
        M12 = u.x*u.y*(1-c)-u.z*s;
        M13 = u.x*u.z*(1-c)+u.y*s;

        M21 = u.x*u.y*(1-c)+u.z*s;
        M22 = u.y*u.y*(1-c)+c;
        M23 = u.y*u.z*(1-c)-u.x*s;

        M31 = u.x*u.z*(1-c)-u.y*s;
        M32 = u.y*u.z*(1-c)+u.x*s;
        M33 = u.z*u.z*(1-c)+c;

//      qx,qy,qz - x,y,z coordinates of the point given by vec after rotation        
        qx = M11*vec.x + M12*vec.y + M13*vec.z + vec1.x;
        qy = M21*vec.x + M22*vec.y + M23*vec.z + vec1.y;
        qz = M31*vec.x + M32*vec.y + M33*vec.z + vec1.z;
        
        vectorResult = new AtomVector(qx, qy, qz);		
		return vectorResult;		
	}
	
	protected StringBuilder makeAllAtomsXYZ(final ArrayList<Molecule> molecules, final HashMap<Integer, MolecularModel> molModelMap, final int atomsTotal){
		
		MolecularModel molModel;
		StringBuilder sb = new StringBuilder();		
				
	     sb.append(atomsTotal)
	       .append("\n\n");
				
		 for (int i = 0; i < molecules.size(); i++) {
					
			molModel = molModelMap.get(molecules.get(i).molModelKey);
					
			for (int j = 0; j < molecules.get(i).numbOfCoulombSites; j++) {
						
				sb.append(molModel.sitesCoulomb.get(j).labelSiteCoulomb).append("\t")
				  .append(molecules.get(i).molCoulombSites.get(j).x).append("\t")
				  .append(molecules.get(i).molCoulombSites.get(j).y).append("\t")
				  .append(molecules.get(i).molCoulombSites.get(j).z)
				  .append("\n");
			}
		}
		
		sb.append("\n");
				
		return sb;
	
	}
	
	protected StringBuilder getStringBuilderAllAtomsXYZ(){
		
		return makeAllAtomsXYZ(molecules, molModelMap, 3*molecules.size());
	}
	
	
	protected void writeXYZ_VMD_File(final String fileName, final ArrayList<Molecule> molecules, final HashMap<Integer, MolecularModel> molModelMap, final int atomsTotal){
				
		MolecularModel molModel;
		
		BufferedWriter bw;
		try {
			bw = new BufferedWriter(new FileWriter("./output_files/" + fileName));
			
			try {
				
				bw.write( atomsTotal + "\n\n");
				
				for (int i = 0; i < molecules.size(); i++) {
					
					molModel = molModelMap.get(molecules.get(i).molModelKey);
					
					for (int j = 0; j < molecules.get(i).numbOfCoulombSites; j++) {
					
						bw.write(String.format("%1s \t %10.5f \t %10.5f \t %10.5f  \n" , 
													molModel.sitesCoulomb.get(j).labelSiteCoulomb,
													molecules.get(i).molCoulombSites.get(j).x,       
													molecules.get(i).molCoulombSites.get(j).y,       
													molecules.get(i).molCoulombSites.get(j).z
								               )
						);			
					}
				}
				
			} catch (IOException e) {
				e.printStackTrace();
			} 
			
			try {
				bw.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
			
		} catch (IOException e1) {
			e1.printStackTrace();
		}		
	}
	

	protected void readXYZ_File(final String fileName, ArrayList<AtomVector> array){
		
		int lineCounter = 0;
		
		try {
			BufferedReader br = new BufferedReader(new FileReader("./input_files/" + fileName));			
			ArrayList lines = new ArrayList();
			
			try {				
				lineCounter = 0;
				
				for(String line = br.readLine();line != null;line = br.readLine()) {
					
					lineCounter++;					
					line = line.trim();					
					String[] fields = line.split("\\s+");
					
					AtomVector vec = new AtomVector(Double.parseDouble(fields[0]), 
							                        Double.parseDouble(fields[1]), 
							                        Double.parseDouble(fields[2])
							                        );  				
					array.add(vec);					
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
			
			try {
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}	
	}
	
	protected void readXYZ_File(final String fileName){
		
		int lineCounter = 0;
		
		try {
			BufferedReader br = new BufferedReader(new FileReader("./input_files/" + fileName));
			
			ArrayList lines = new ArrayList();
			
			try {
				
				lineCounter = 0;
				
				for(String line = br.readLine();line != null;line = br.readLine()) {
					
					lineCounter++;					
					line = line.trim();					
					String[] fields = line.split("\\s+");
				
					System.out.print(lineCounter + "\t"); 
		            System.out.print("fields[0] = " + fields[0] + "\t"); 
		            System.out.print("fields[1] = " + fields[1] + "\t"); 
		            System.out.println("fields[2] = " + fields[2]); 
					
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}		
	}
	
	
	protected Integer nextKeyInMap(final HashMap<Integer, MolecularModel> map){
		
		Integer keyMap;

		if (map.isEmpty()) {
			keyMap = 0;
		} else {
			keyMap = Collections.max(molModelMap.keySet()) + 1;
		}
		return keyMap;
	}
	
	//  Large promoter molecule THF model setup:
	protected MolecularModel setupTHFModel(){
		
		ArrayList<LJSite> sitesLJ;
		AtomVector ljSiteVector;
		LJSite ljSite;
		
		sitesLJ = new ArrayList<LJSite>();
		
		ljSiteVector = new AtomVector(0.000012,	-1.192456, -0.314052);		
		ljSite = new LJSite("O", 0.1700, 3.000012, ljSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesLJ.add(ljSite);
		
		ljSiteVector = new AtomVector(1.120439,	-0.465886, 0.167894);		
		ljSite = new LJSite("C", 0.1094, 3.399670, ljSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesLJ.add(ljSite);
		
		ljSiteVector = new AtomVector(-1.120429, -0.465909,	0.167891);		
		ljSite = new LJSite("C", 0.1094, 3.399670, ljSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesLJ.add(ljSite);
		
		ljSiteVector = new AtomVector(0.771947,	1.008301, -0.053734);		
		ljSite = new LJSite("C", 0.1094, 3.399670, ljSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesLJ.add(ljSite);
		
		ljSiteVector = new AtomVector(-0.771968, 1.008287, -0.053731);		
		ljSite = new LJSite("C", 0.1094, 3.399670, ljSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesLJ.add(ljSite);
		
		ljSiteVector = new AtomVector(1.997366,	-0.803244, -0.3763);		
		ljSite = new LJSite("H", 0.0157, 2.649533, ljSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesLJ.add(ljSite);
		
		ljSiteVector = new AtomVector(-1.99735,	-0.803282, -0.376305);		
		ljSite = new LJSite("H", 0.0157, 2.649533, ljSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesLJ.add(ljSite);
		
		ljSiteVector = new AtomVector(1.259694,	-0.672439, 1.234558);		
		ljSite = new LJSite("H", 0.0157, 2.649533, ljSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesLJ.add(ljSite);
		
		ljSiteVector = new AtomVector(-1.259685, -0.672467,	1.234554);		
		ljSite = new LJSite("H", 0.0157, 2.649533, ljSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesLJ.add(ljSite);
		
		ljSiteVector = new AtomVector(1.155732,	1.349436, -1.011001);		
		ljSite = new LJSite("H", 0.0157, 2.649533, ljSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesLJ.add(ljSite);
		
		ljSiteVector = new AtomVector(-1.155763, 1.34942, -1.010994);		
		ljSite = new LJSite("H", 0.0157, 2.649533, ljSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesLJ.add(ljSite);
		
		ljSiteVector = new AtomVector(1.190595,	1.641746, 0.723987);		
		ljSite = new LJSite("H", 0.0157, 2.649533, ljSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesLJ.add(ljSite);
		
		ljSiteVector = new AtomVector(-1.190624, 1.64172, 0.723997);		
		ljSite = new LJSite("H", 0.0157, 2.649533, ljSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesLJ.add(ljSite);
		
		
		ArrayList<CoulombSite> sitesCoulomb;
		AtomVector coulSiteVector;	
		CoulombSite coulSite;
		
		sitesCoulomb = new ArrayList<CoulombSite>();
		
		coulSiteVector = new AtomVector(0.000012,	-1.192456, -0.314052);		
		coulSite = new CoulombSite("O", -0.42513, coulSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesCoulomb.add(coulSite);
		
		coulSiteVector = new AtomVector(1.120439,	-0.465886, 0.167894);		
		coulSite = new CoulombSite("C", 0.171469, coulSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesCoulomb.add(coulSite);
		
		coulSiteVector = new AtomVector(-1.120429, -0.465909,	0.167891);		
		coulSite = new CoulombSite("C", 0.171478, coulSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesCoulomb.add(coulSite);
		
		coulSiteVector = new AtomVector(0.771947,	1.008301, -0.053734);		
		coulSite = new CoulombSite("C", 0.024169, coulSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesCoulomb.add(coulSite);
		
		coulSiteVector = new AtomVector(-0.771968, 1.008287, -0.053731);		
		coulSite = new CoulombSite("C", 0.024154, coulSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesCoulomb.add(coulSite);
		
		coulSiteVector = new AtomVector(1.997366,	-0.803244, -0.3763);		
		coulSite = new CoulombSite("H", 0.022566, coulSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesCoulomb.add(coulSite);
		
		coulSiteVector = new AtomVector(-1.99735,	-0.803282, -0.376305);		
		coulSite = new CoulombSite("H", 0.022564, coulSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesCoulomb.add(coulSite);
		
		coulSiteVector = new AtomVector(1.259694,	-0.672439, 1.234558);		
		coulSite = new CoulombSite("H", 0.007663, coulSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesCoulomb.add(coulSite);
		
		coulSiteVector = new AtomVector(-1.259685, -0.672467,	1.234554);		
		coulSite = new CoulombSite("H", 0.007661, coulSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesCoulomb.add(coulSite);
		
		coulSiteVector = new AtomVector(1.155732,	1.349436, -1.011001);		
		coulSite = new CoulombSite("H", 0.005161, coulSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesCoulomb.add(coulSite);
		
		coulSiteVector = new AtomVector(-1.155763, 1.34942, -1.010994);		
		coulSite = new CoulombSite("H", 0.005164, coulSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesCoulomb.add(coulSite);
		
		coulSiteVector = new AtomVector(1.190595,	1.641746, 0.723987);		
		coulSite = new CoulombSite("H", -0.018462, coulSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesCoulomb.add(coulSite);
		
		coulSiteVector = new AtomVector(-1.190624, 1.64172, 0.723997);		
		coulSite = new CoulombSite("H", -0.018458, coulSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesCoulomb.add(coulSite);
		
		
		MolecularModelType molModelTypeProm_L_THF = MolecularModelType.PROMOTOR_LARGE;		
		MolecularModel promModelProm_L_THF = new MolecularModel("THF", molModelTypeProm_L_THF, sitesLJ.size(), sitesLJ, sitesCoulomb.size(), sitesCoulomb);
		
		return promModelProm_L_THF;		
	}

	//  Large promoter molecule Acetylene model setup:
	protected MolecularModel setupAcetyleneModel(){

		ArrayList<LJSite> sitesLJ = new ArrayList<LJSite>();	
		
		AtomVector ljSiteVector = new AtomVector(0.000000, 0.000000, 0.605696);		
		LJSite ljSite = new LJSite("C", 0.08600, 3.39967, ljSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesLJ.add(ljSite);
		
		ljSiteVector = new AtomVector(0.000000,	0.000000, -0.605696);		
		ljSite = new LJSite("C", 0.08600, 3.39967, ljSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesLJ.add(ljSite);
		
		ljSiteVector = new AtomVector(0.000000,	0.000000, 1.667240);		
		ljSite = new LJSite("H", 0.01570, 2.64953, ljSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesLJ.add(ljSite);
		
		ljSiteVector = new AtomVector(0.000000, 0.000000, -1.667240);		
		ljSite = new LJSite("H", 0.01570, 2.64953, ljSiteVector);		// epsilon in kcal/mol, sigma in angstroms 	
		sitesLJ.add(ljSite);
		
		
		ArrayList<CoulombSite> sitesCoulomb = new ArrayList<CoulombSite>();
		
		AtomVector coulSiteVector = new AtomVector(0.000000, 0.000000, 0.605696);
		CoulombSite coulSite = new CoulombSite("C", -0.263874, coulSiteVector);
		sitesCoulomb.add(coulSite);
		
		coulSiteVector = new AtomVector(0.000000,	0.000000, -0.605696);
		coulSite = new CoulombSite("C", -0.263874, coulSiteVector);
		sitesCoulomb.add(coulSite);
		
		coulSiteVector = new AtomVector(0.000000,	0.000000, 1.667240);
		coulSite = new CoulombSite("H", 0.263874, coulSiteVector);
		sitesCoulomb.add(coulSite);
		
		coulSiteVector = new AtomVector(0.000000, 0.000000, -1.667240);
		coulSite = new CoulombSite("H", 0.263874, coulSiteVector);
		sitesCoulomb.add(coulSite);
		
		MolecularModelType molModelTypeProm_L_ACT = MolecularModelType.PROMOTOR_LARGE;		
		MolecularModel promModelProm_L_ACT = new MolecularModel("Acetylene", molModelTypeProm_L_ACT, sitesLJ.size(), sitesLJ, sitesCoulomb.size(), sitesCoulomb);
		
		return promModelProm_L_ACT;
	}

//	Hydrogen molecule model setup
	protected MolecularModel setupHydrogenModel(){
//		Saman Alavi; J A Ripmeester; D D Klug
//		Molecular dynamics simulations of binary structure H hydrogen and methyl-tert-butylether clathrate hydrates.
//		The Journal of chemical physics 2006;124(20):204707
//		http://link.aip.org/link/doi/10.1063/1.2199850

		ArrayList<LJSite> sitesLJ = new ArrayList<LJSite>();
		
		AtomVector ljSiteVector = new AtomVector(0.0, 0.0, 0.0);	
		
		LJSite ljSite = new LJSite("M", 0.06811884971816, 3.038, ljSiteVector);		// epsilon in kcal/mol, sigma in angstroms 		
		sitesLJ.add(ljSite);
		
		ArrayList<CoulombSite> sitesCoulomb = new ArrayList<CoulombSite>();
		
		AtomVector coulSiteVector = new AtomVector(0.0, 0.0, 0.0);
		CoulombSite coulSite = new CoulombSite("M", -0.9864, coulSiteVector);
		sitesCoulomb.add(coulSite);
		
		coulSiteVector = new AtomVector(0.7414/2.0, 0.0, 0.0);
		coulSite = new CoulombSite("H", 0.4932, coulSiteVector);
		sitesCoulomb.add(coulSite);
		
		coulSiteVector = new AtomVector(-0.7414/2.0, 0.0, 0.0);
		coulSite = new CoulombSite("H", 0.4932, coulSiteVector);
		sitesCoulomb.add(coulSite);
		
		MolecularModelType molModelTypeProm_SM_H2 = MolecularModelType.PROMOTOR_SMALL_MEDIUM;		
		MolecularModel hydrogenModelProm_SM_H2 = new MolecularModel("Hydrogen", molModelTypeProm_SM_H2, sitesLJ.size(), sitesLJ, sitesCoulomb.size(), sitesCoulomb);
		
		return hydrogenModelProm_SM_H2;
	}

//		Setup of TIP3P model of water
	protected MolecularModel setupWaterModelTIP3P(final double alw, final double dHOw){
		
        // al is an angle in radians between OH bonds in water molecule
        double al = (alw/2.0)*Math.PI/180.0; // bisector angle O-H-O in radians      

		ArrayList<LJSite> sitesLJ = new ArrayList<LJSite>();
        
		AtomVector ljSiteVector = new AtomVector(0.0, 0.0, 0.0);		
		LJSite ljSite = new LJSite("O", 0.1521033, 3.15061, ljSiteVector);		// epsilon in kcal/mol, sigma in angstroms 
		sitesLJ.add(ljSite);
		
		final double H1x = dHOw*Math.sin(al);
		final double H1y = dHOw*Math.cos(al);
		final double H1z = 0.0;

		final double H2x = -H1x;
		final double H2y =  H1y;
		final double H2z =  0.0;
		
		
		ArrayList<CoulombSite> sitesCoulomb = new ArrayList<CoulombSite>();	
		
		AtomVector coulSiteVector = new AtomVector(0.0, 0.0, 0.0);
		CoulombSite coulSite = new CoulombSite("O", -0.8340, coulSiteVector);
		sitesCoulomb.add(coulSite);
	
		coulSiteVector = new AtomVector(H1x, H1y, H1z);
		coulSite = new CoulombSite("H", 0.4170, coulSiteVector);
		sitesCoulomb.add(coulSite);
		
		coulSiteVector = new AtomVector(H2x, H2y, H2z);
		coulSite = new CoulombSite("H", 0.4170, coulSiteVector);
		sitesCoulomb.add(coulSite);	
		
		MolecularModelType molModelTypeWater = MolecularModelType.WATER;		
		MolecularModel waterModelWater = new MolecularModel("TIP3P", molModelTypeWater, sitesLJ.size(), sitesLJ, sitesCoulomb.size(), sitesCoulomb);
//		System.out.println(waterModelWater);
		
		return waterModelWater;		
	}
	
//	Setup of TIP3P model of water
protected MolecularModel setupWaterModelSPCE(final double alw, final double dHOw){
	
    // al is an angle in radians between OH bonds in water molecule
    double al = (alw/2.0)*Math.PI/180.0; // bisector angle O-H-O in radians      

	ArrayList<LJSite> sitesLJ = new ArrayList<LJSite>();
    
	AtomVector ljSiteVector = new AtomVector(0.0, 0.0, 0.0);		
	LJSite ljSite = new LJSite("O", 0.155354, 3.166, ljSiteVector);		// epsilon in kcal/mol, sigma in angstroms 
	sitesLJ.add(ljSite);
	
	final double H1x = dHOw*Math.sin(al);
	final double H1y = dHOw*Math.cos(al);
	final double H1z = 0.0;

	final double H2x = -H1x;
	final double H2y =  H1y;
	final double H2z =  0.0;
	
	
	ArrayList<CoulombSite> sitesCoulomb = new ArrayList<CoulombSite>();	
	
	AtomVector coulSiteVector = new AtomVector(0.0, 0.0, 0.0);
	CoulombSite coulSite = new CoulombSite("O", -0.8476, coulSiteVector);
	sitesCoulomb.add(coulSite);

	coulSiteVector = new AtomVector(H1x, H1y, H1z);
	coulSite = new CoulombSite("H", +0.4238, coulSiteVector);
	sitesCoulomb.add(coulSite);
	
	coulSiteVector = new AtomVector(H2x, H2y, H2z);
	coulSite = new CoulombSite("H", +0.4238, coulSiteVector);
	sitesCoulomb.add(coulSite);	
	
	MolecularModelType molModelTypeWater = MolecularModelType.WATER;		
	MolecularModel waterModelWater = new MolecularModel("SPC/E", molModelTypeWater, sitesLJ.size(), sitesLJ, sitesCoulomb.size(), sitesCoulomb);
//	System.out.println(waterModelWater);
	
	return waterModelWater;		
}

protected void writeRDF_File(final String fileName, final AtomVector box, final int numbWaterMols, final int numbMCSimulations, int [] g){
	
	double vol, r, rdf;
	
	vol = box.x * box.y * box.z;
	
	BufferedWriter bw;
	try {
		bw = new BufferedWriter(new FileWriter("./output_files/" + fileName));
		
		try {
			
			r = 0.0;
			
			for (int i = 0; i < g.length; i++) {
				
				rdf =  vol * g[i] / (4.0 * Math.PI * Math.pow((i*dR - 0.5*dR), 2) * dR * numbWaterMols * numbMCSimulations );
				r = r + dR;

				if ( r <= 9.6)
					bw.write(String.format("%10.5f \t %10.5f \n" , r, rdf ) );
			}
			
		} catch (IOException e) {
			e.printStackTrace();
		} 
		
		try {
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}		
		
	} catch (IOException e1) {
		e1.printStackTrace();
	}		
}

	protected void normalizeRDF(final AtomVector box, final int numbWaterMols,
			final int numbMCSimulations, int[] g) {

		double vol, r, rdf;
		Point2D pt = null;
		
		gr.clear();

		vol = box.x * box.y * box.z;
		r = 0.0;

		for (int i = 0; i < g.length; i++) {

			rdf = vol * g[i]
					/ (4.0 * Math.PI * Math.pow((i * dR - 0.5 * dR), 2) * dR * numbWaterMols * numbMCSimulations);
			r = r + dR;

			if (r <= 9.6) {
				pt = new Point2D.Double(r, rdf);
				gr.add(pt);
			}

		}
	}		

	
	protected final void calculateRDF(final AtomVector box,
									  final ArrayList<Molecule> molecules, 
			                          final int numbWaterMols,
									  int [] g									  
			                         ){

		int jatom, index;
		double r, rSq;
		
		AtomVector atomI = new AtomVector();
		AtomVector atomJ = new AtomVector();
		AtomVector vecIJ = new AtomVector();
		AtomVector pbc   = new AtomVector();
		
		LinkedList<Double> listBox = new LinkedList<Double>();
		listBox.add(box.x);
		listBox.add(box.y);
		listBox.add(box.z);
		double halfBox = 0.5*Collections.min(listBox);
		
		Random random = new Random();
		
		int start = 0;
	    int end = numbWaterMols - 1;		    
	    int iMol = RandomInteger(start, end, random);   // generate random index of molecule
	    
	    atomI =  molecules.get(iMol).molCoulombSites.get(0);
	    
	    for (int jMol=0; jMol < numbWaterMols; jMol++){
	    	
	    	if (iMol != jMol){
	    	
				start = 1;
				end = 2;
				jatom = RandomInteger(start, end, random);

				atomJ = molecules.get(jMol).molCoulombSites.get(jatom);

				vecIJ = atomI.AtomVectorSubtract(atomJ);

				pbc.x = box.x * Math.rint(vecIJ.x / box.x);
				pbc.y = box.y * Math.rint(vecIJ.y / box.y);
				pbc.z = box.z * Math.rint(vecIJ.z / box.z);

				vecIJ = vecIJ.AtomVectorSubtract(pbc);
				rSq = vecIJ.AtomVectorLengthSquared();

				if (rSq < halfBox * halfBox) {

					r = Math.sqrt(rSq);
					index = (int) Math.rint(r / dR);
					g[index] += 1;
				}
					
	    		
	    	}	    	
	    }
		
	}
	
	protected void simRun(){
		
		int iMol, iRan;
		int start, end;
		Random random = new Random();
//		random.setSeed(123456789);
		
		StringBuilder sb = new StringBuilder();
		
		for (int i=0; i < numbMCSimulations; i++){
			
			if ( !go ) break;
			
			stepMC = i;
			
			start = 0;
		    end = numbOfMolsActual - 1;		    
		    iMol = RandomInteger(start, end, random);   // generate random index of molecule
		    
			start = 1;
		    end = 2*numbOfMolsActual + 3;		    
		    iRan = RandomInteger(start, end, random);
			
			engineMonteCarlo(i, iMol, iRan);	
			
			calculateRDF(box, molecules, numbWaterMols, g);
			
			if ((stepMC % (freqOfSysOut-5)) == 0) {
				simulOutputData.setData(i, enTot, enTot, box, gr, sb);
			}
			
		}
		
		writeRDF_File("output_simulation_g_Ow_Hw.xyz", box, numbWaterMols, numbMCSimulations, g);
		writeXYZ_VMD_File("output_simulation_water_framework_VMD.xyz", molecules, molModelMap, 3*molecules.size());
		
	}
	
	protected void rdf(){
		calculateRDF(box, molecules, numbWaterMols, g);
	}

}
