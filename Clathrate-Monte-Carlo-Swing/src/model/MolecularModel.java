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

public class MolecularModel {
	
	public String molModelName;	
	public MolecularModelType molModelType;	
	public int numbOfLJSites;
	public int numbOfCoulombSites;
	
	public ArrayList<LJSite> sitesLJ;
	public ArrayList<CoulombSite> sitesCoulomb;
	
	MolecularModel(  String molModelName,
					 MolecularModelType molModelType,
			         int numbOfLJSites,
			         ArrayList<LJSite> sitesLJ,
	                 int numbOfCoulombSites,	                 
	                 ArrayList<CoulombSite> sitesCoulomb			
				  ){
		
		this.molModelName = molModelName;
		this.molModelType = molModelType;
		this.numbOfLJSites = numbOfLJSites;
		this.numbOfCoulombSites = numbOfCoulombSites;
		this.sitesLJ = sitesLJ;
		this.sitesCoulomb = sitesCoulomb;
	}

	@Override
	public String toString() {
		return "MolecularModel [molModelName=" + molModelName + ", molModelType=" + molModelType + ", numbOfLJSites="
				+ numbOfLJSites + ", sitesLJ=" + sitesLJ + ", numbOfCoulombSites=" + numbOfCoulombSites 
				+ ", sitesCoulomb=" + sitesCoulomb + "]";
	}


}
