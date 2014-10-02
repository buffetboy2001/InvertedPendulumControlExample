/**
 * 
 */
package com.github.buffetboy2001;

import javax.swing.JFrame;

import org.mitre.caasd.jlcl.components.FixedStepBackwardDiffentiatorArguments;
import org.mitre.caasd.jlcl.components.FixedStepBackwardDifferentiator;
import org.mitre.caasd.jlcl.components.FixedStepEulerIntegration;
import org.mitre.caasd.jlcl.components.FixedStepIntegrationArguments;
import org.mitre.caasd.jlcl.components.FixedStepIntegrator;
import org.mitre.caasd.jlcl.components.PID;
import org.mitre.caasd.jlcl.components.ProportionalIntegralDerivativeArguments;
import org.mitre.caasd.jlcl.interfaces.IPID;
import org.mitre.caasd.jlcl.interfaces.IPidEvaluationArguments;
import org.opensourcephysics.frames.PlotFrame;
import org.opensourcephysics.numerics.Euler;
import org.opensourcephysics.numerics.ODE;
import org.opensourcephysics.numerics.ODESolver;

/**
 * @author SBOWMAN
 *
 */
public class InvertedPendulumODE implements ODE {

	private double time = 0; // actually only needed for visualization
	private double[] propogationState = new double[4], previousState = new double[4], previousPreviousState = new double[4]; // x, xdot, phi, phidot
	private ODESolver solver = new Euler(this);
	private boolean isAtPhiMax = false;
	private IPID<Double,FixedStepIntegrator<Double,FixedStepIntegrationArguments<Double>>,FixedStepIntegrationArguments<Double>,FixedStepBackwardDifferentiator<Double>,FixedStepBackwardDiffentiatorArguments<Double>> pidController = null;
	private FixedStepIntegrationArguments<Double> integrationArgs;
	private FixedStepBackwardDiffentiatorArguments<Double> differentiatorArgs;
	private double controlForce = 0;
	
	// Constants of the solution
	final static double stepSize = 0.01; // seconds
	final double propGain = 10;
	final double intGain = 0.5;
	final double derivGain = 10;
	final double cartMass = 0.5; // kg
	final double pendulumMass = 0.2; // kg
	final double coefOfFriction = 0.1; // N/meter/sec
	final double lengthOfPendulum = 0.3; // meter
	final double momentOfIntertiaPendulum = 0.006; // kg m^2
	final double gravityConst = 9.81; // meters/sec^2
	final double phiMax = Math.PI/2;
	
	// Derived constants
	final double lengthOfPendulumSquared = lengthOfPendulum*lengthOfPendulum;
	
	// Visualization
	private PlotFrame frame;
	
	public InvertedPendulumODE() {
		
		// Hardcode initial state
		propogationState[0] = 0;
		propogationState[1] = 0;
		propogationState[2] = 1/57.3; // radians, positive is left of center
		propogationState[3] = 0;
		previousState = propogationState;
		previousPreviousState = previousPreviousState;
		
		// Set solver step size (time, seconds)
		solver.initialize(stepSize);
		
		// Instantiate controller
		integrationArgs = new FixedStepIntegrationArguments<Double>(Double.class,solver.getStepSize());
		FixedStepIntegrator<Double,FixedStepIntegrationArguments<Double>> integrator = new FixedStepEulerIntegration<Double>(Double.class,integrationArgs);
		differentiatorArgs = new FixedStepBackwardDiffentiatorArguments<Double>(Double.class,solver.getStepSize());
		FixedStepBackwardDifferentiator<Double> differentiator = new FixedStepBackwardDifferentiator<Double>(Double.class);
		pidController = new PID<Double,FixedStepIntegrator<Double,FixedStepIntegrationArguments<Double>>,FixedStepIntegrationArguments<Double>,FixedStepBackwardDifferentiator<Double>,FixedStepBackwardDiffentiatorArguments<Double>>(Double.class,propGain,integrator,intGain,differentiator,derivGain);
		
		// Prepare visualization
		createPlotFrame();
	}
	
	/* (non-Javadoc)
	 * @see org.opensourcephysics.numerics.ODE#getState()
	 */
	public double[] getState() {
		return this.propogationState;
	}

	/* (non-Javadoc)
	 * @see org.opensourcephysics.numerics.ODE#getRate(double[], double[])
	 */
	public void getRate(double[] state, double[] rate) {
		// Force a limit on phi
		if (!isAtPhiMax && Math.abs(state[2]) > phiMax) {
			isAtPhiMax = true;
			state[2] = Math.signum(state[2])*phiMax; // limit phi to maximum
			state[3] = 0; // by definition no longer rotating
		}
				
		// debug
		outputToConsole("state",state); 
		
		// Rate Calculation based on the current state
		double denominator = momentOfIntertiaPendulum*(cartMass+pendulumMass)+(pendulumMass*cartMass*lengthOfPendulumSquared);
		rate[0] = state[1]; // velocity of cart
		rate[1] = state[1] * (-coefOfFriction*(momentOfIntertiaPendulum+pendulumMass*lengthOfPendulumSquared)/denominator) + state[2] * (pendulumMass*pendulumMass*gravityConst*lengthOfPendulumSquared/denominator);
		rate[2] = state[3]; // rotation velocity of pendulum
		rate[3] = state[1] * (-pendulumMass*lengthOfPendulum*coefOfFriction/denominator) + state[2] * (pendulumMass*gravityConst*lengthOfPendulum*(cartMass+pendulumMass)/denominator);
		if (isAtPhiMax) {
			rate[3] = 0; // pendulum is no longer accelerating
		}
		
		// Add control mechanism and control force effects
		rate[1] += controlForce*(momentOfIntertiaPendulum+pendulumMass*lengthOfPendulumSquared)/denominator;
		rate[3] += controlForce*(pendulumMass*lengthOfPendulum)/denominator;
		
		// Hold on to the latest state for the feedback loop
		previousPreviousState = previousState;
		previousState = propogationState;
		propogationState = state;

		// debug
//		outputToConsole("rate",rate);

	}
	
	public void doSolverStep() {
		// Calculate the pendulum equations of motion
		time += stepSize; // only for visualization
		solver.step();

		// Calculate the controller output for the current state. Feedback phi only and the reference signal is 0. So, the error value is -phi.
		IPidEvaluationArguments<Double,FixedStepIntegrationArguments<Double>,FixedStepBackwardDiffentiatorArguments<Double>> pidArgs = new ProportionalIntegralDerivativeArguments<Double,FixedStepIntegrationArguments<Double>,FixedStepBackwardDiffentiatorArguments<Double>>(Double.class);
		differentiatorArgs.updatePreviousFunctionValue(-previousState[2]);
		pidArgs.setIntegratorEvaluationArgs(integrationArgs);
		pidArgs.setDifferentiatorEvaluationArgs(differentiatorArgs);
		pidArgs.updateErrorSignalValue(-propogationState[2]);
		controlForce = pidController.evaluate(pidArgs);
		
		// update plot
		updatePlot(time, propogationState[2]*180/Math.PI);
	}

	private void createPlotFrame() {
	    frame = new PlotFrame("Time (sec)", "Phi", "Controller Behavior"); 
	    frame.setSize(400, 400);
	    frame.setVisible(true);
	    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setAutoscaleX(true);
		frame.setAutoscaleY(true);
	}
	
	private void updatePlot(double x, double y) {
		frame.append(0, x, y);
	}

	public static void main(String[] args) {
		final double finalStopTime = 4; // seconds
		InvertedPendulumODE ip = new InvertedPendulumODE(); 
		for (int i = 0; i<finalStopTime/ip.stepSize;i++)
			ip.doSolverStep();
	}
	
	public static void outputToConsole(String prefix, double[] array) {
		System.out.println(prefix + ":\n\tx-pos: " + array[0] + "\n\tx-dot: " + array[1] + "\n\tphi: " + array[2]*57.3 + "\n\tphi-dot: " + array[3]);		
	}
}
