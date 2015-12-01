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

import java.util.ArrayList;

public class Molecule {
	
	public Integer molModelKey;	
	
	public int numbOfLJSites;
	public int numbOfCoulombSites;
	
	public AtomVector molCenterOfMass;
	public MolEulerAngles molEulerAngles;
	
	public ArrayList<AtomVector> molLJSites;
	public ArrayList<AtomVector> molCoulombSites;
	
	Molecule(	Integer molModelKey,
			 	int numbOfLJSites,
			 	int numbOfCoulombSites,
			 	AtomVector molCenterOfMass,
			 	MolEulerAngles molEulerAngles,
			 	ArrayList<AtomVector> molLJSites,
			 	ArrayList<AtomVector> molCoulombSites
			 ){
		
		this.molModelKey = molModelKey;
		this.numbOfLJSites = numbOfLJSites;
		this.numbOfCoulombSites = numbOfCoulombSites;
		this.molCenterOfMass = molCenterOfMass;
		this.molEulerAngles = molEulerAngles;
		this.molLJSites = molLJSites;
		this.molCoulombSites = molCoulombSites;
	}

	@Override
	public String toString() {
		return "Molecule [molModelKey=" + molModelKey + ", numbOfLJSites=" + numbOfLJSites + ", numbOfCoulombSites="
				+ numbOfCoulombSites + ", molCenterOfMass=" + molCenterOfMass + ", molEulerAngles=" + molEulerAngles
				+ ", molLJSites=" + molLJSites + ", molCoulombSites=" + molCoulombSites + "]";
	}

}


