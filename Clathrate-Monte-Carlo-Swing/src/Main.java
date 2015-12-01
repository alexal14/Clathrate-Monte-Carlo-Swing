import java.awt.EventQueue;

import model.*;
import view.*;
import controller.*;

public class Main {

	public static void main(String[] args) {
		
		try {
			Model model = new Model();
			View view = new View();
			Controller controller = new Controller(model, view);
			view.setVisible(true);
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}

}





