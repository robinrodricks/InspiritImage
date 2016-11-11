package ru.inspirit.steering.behavior {
	import ru.inspirit.steering.SteerVector3D;
	import ru.inspirit.steering.Vehicle;
	import ru.inspirit.steering.behavior.AbstractBehavior;

	/**
	 * @author Eugene Zatepyakin
	 */
	public class Flee extends AbstractBehavior 
	{
		public var target:SteerVector3D = new SteerVector3D();
		
		public function Flee(target:SteerVector3D) 
		{
			this.target = target;
		}
		
		override public function apply(veh:Vehicle):void
		{
			accumulator.setDiff(veh.position, target);
			
			var goalLength:Number = 1.1 * veh.velocity.approximateLength();
			accumulator.setApproximateTruncate(goalLength);
			accumulator.setDiff(accumulator, veh.velocity);
			accumulator.setApproximateTruncate(veh.maxForce);
			
			veh.applyGlobalForce(accumulator);
		}
	}
}
