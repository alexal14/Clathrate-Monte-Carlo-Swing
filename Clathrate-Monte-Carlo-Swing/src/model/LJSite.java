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

public class LJSite {
	
	public String labelLJSite;	
	public double epsilon;
	public double sigma;	
	public AtomVector siteLJ;
	
	LJSite(String labelLJSite, double epsilon, double sigma, AtomVector siteLJ){
		
		this.labelLJSite = labelLJSite;
		this.epsilon = epsilon;
		this.sigma = sigma;
		this.siteLJ = siteLJ;
	}

	@Override
	public String toString() {
		return "LJSite [labelLJSite=" + labelLJSite + ", epsilon=" + epsilon + ", sigma=" + sigma + ", siteLJ=" + siteLJ + "]";
	}
	
	

}
