package ru.inspirit.steering.behavior 
{
	import ru.inspirit.steering.SteerVector3D;
	import ru.inspirit.steering.Vehicle;

	/**
	 * @author Eugene Zatepyakin
	 */
	public class AbstractBehavior 
	{		
		public var accumulator:SteerVector3D = new SteerVector3D();
		
		public var next:AbstractBehavior;
		
		public function apply(veh:Vehicle):void {}
		
		public function nearby(veh:Vehicle, other:Vehicle, distance:Number, angleCosine:Number):Number
		{			
			var dx:Number = other.position.x - veh.position.x;
			var dy:Number = other.position.y - veh.position.y;
			var dz:Number = other.position.z - veh.position.z;
			
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
			
			var approximateDistance:Number = a * 0.9375 + (b + c) * 0.375;
			
			if(approximateDistance > distance)
			{
				return 0;
			}
			
			var dot:Number = dx * veh.forward.x + dy * veh.forward.y + dz * veh.forward.z;
			
			return dot >= angleCosine ? approximateDistance : 0;
		}
	}
}
