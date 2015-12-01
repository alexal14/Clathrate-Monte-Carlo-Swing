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

public class MolEulerAngles {
	
	public double alpha;  // or fi, defined in the range [-pi,+pi]
	public double beta;   // or theta, defined in the range [0,+pi]
	public double gamma;  // or psi, defined in the range [-pi,+pi] 
	
	MolEulerAngles(){
		this.alpha = 0.0;
		this.beta = 0.0;
		this.gamma = 0.0;
	}
	
	MolEulerAngles(double alpha, double beta, double gamma){
		this.alpha = alpha;
		this.beta = beta;
		this.gamma = gamma;
	}

	@Override
	public String toString() {
		return "MolEulerAngles [alpha=" + alpha + ", beta=" + beta + ", gamma=" + gamma + "]";
	}

}
