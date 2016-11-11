package ru.inspirit.steering.behavior 
{
	import ru.inspirit.steering.SteerVector3D;
	import ru.inspirit.steering.Vehicle;
	import ru.inspirit.steering.behavior.AbstractBehavior;

	/**
	 * @author Eugene Zatepyakin
	 */
	public final class Arrival extends AbstractBehavior 
	{
		public var target:SteerVector3D;
		public var slowingDistance:Number = 7;
		
		public function Arrival(target:SteerVector3D, slowingDistance:Number = 7)
		{
			this.target = target;
			this.slowingDistance = slowingDistance;
		}
		
		override public function apply(veh:Vehicle):void
		{
			var dx:Number = target.x - veh.position.x;
			var dy:Number = target.y - veh.position.y;
			var dz:Number = target.z - veh.position.z;
			
			var a:Number = dx;
			if (a < 0) a = -a;
			
			var b:Number = dy;
			if (b < 0) b = -b;
			
			var c:Number = dz;
			if (c < 0) c = -c;
			
			var t:Number;
			if(a < b) {
				t = a;
				a = b;
				b = t;
			}
			if(a < c) {
				t = a;
				a = c;
				c = t;
			}
			
			var distance:Number =  a * 0.9375 + (b + c) * 0.375;
			
			var rampedSpeed:Number = veh.maxSpeed * (distance / slowingDistance);
			var clippedSpeed:Number = (rampedSpeed > veh.maxSpeed ? veh.maxSpeed : rampedSpeed) / distance;
			
			accumulator.x = dx * clippedSpeed - veh.velocity.x;
			accumulator.y = dy * clippedSpeed - veh.velocity.y;
			accumulator.z = dz * clippedSpeed - veh.velocity.z;
			
			veh.applyGlobalForce(accumulator);
		}
	}
}
