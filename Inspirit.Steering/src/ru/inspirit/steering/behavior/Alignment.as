package ru.inspirit.steering.behavior 
{
	import ru.inspirit.steering.Vehicle;
	import ru.inspirit.steering.behavior.AbstractBehavior;

	/**
	 * @author Eugene Zatepyakin
	 */
	public final class Alignment extends AbstractBehavior 
	{
		public var alignDist:Number;
		public var alignAngleCos:Number;
		public var alignStrength:Number;
		
		public var alignList:Vehicle;
		
		public var deltaHeading:Boolean = false;
		
		public function Alignment(alignList:Vehicle = null, alignDist:Number = 9, alignAngle:Number = 90, alignStrength:Number = 0.13, deltaHeading:Boolean = false)
		{
			this.alignList = alignList;
			this.alignDist = alignDist;
			this.alignAngleCos = Math.cos(((alignAngle / 360) * 2.0) * Math.PI);
			this.alignStrength = alignStrength;
			this.deltaHeading = deltaHeading;
		}
		
		override public function apply(veh:Vehicle):void
		{
			var ax:Number = 0, ay:Number = 0, az:Number = 0;
			
			var count:int = 0;
			
			var other:Vehicle = alignList;
			
			while(other) 
			{
				if(other == veh)
				{
					other = other.next;
					continue;
				}
				
				if( nearby(veh, other, alignDist, alignAngleCos) )
				{
					count++;
					
					ax += other.forward.x;
					ay += other.forward.y;
					az += other.forward.z;
				}
				other = other.next;
			}
			
			if(count)
			{
				var f:Number = 1 / count;
				ax *= f;
				ay *= f;
				az *= f;
				
				if (deltaHeading) 
				{
					accumulator.x = ax - veh.forward.x;
					accumulator.y = ay - veh.forward.y;
					accumulator.z = az - veh.forward.z;
				} else {
					accumulator.x = ax;
					accumulator.y = ay;
					accumulator.z = az;
				}
				
				accumulator.setApproximateNormalize();
			
				accumulator.x *= alignStrength;
				accumulator.y *= alignStrength;
				accumulator.z *= alignStrength;
				
				veh.applyGlobalForce(accumulator);
			}
		}
	}
}
