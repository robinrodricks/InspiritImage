package ru.inspirit.steering.behavior 
{
	import ru.inspirit.steering.Vehicle;

	/**
	 * @author Eugene Zatepyakin
	 */
	public final class Separation extends AbstractBehavior 
	{
		public var separateDist:Number;
		public var separateAngleCos:Number;
		public var separateStrength:Number;
		
		public var separateList:Vehicle;
		
		public function Separation(separateList:Vehicle = null, separateDist:Number = 9, separateAngle:Number = 135, separateStrength:Number = 0.95)
		{
			this.separateList = separateList;
			this.separateDist = separateDist;
			this.separateAngleCos = Math.cos(((separateAngle / 360) * 2.0) * Math.PI);
			this.separateStrength = separateStrength;
		}

		override public function apply(veh:Vehicle):void
		{
			var ax:Number = 0, ay:Number = 0, az:Number = 0;
			
			var other:Vehicle = separateList;
			var approximateDistance:Number;
			
			var px:Number = veh.position.x;
			var py:Number = veh.position.y;
			var pz:Number = veh.position.z;
			var dx:Number, dy:Number, dz:Number;
			var f:Number;
			
			while(other) 
			{
				if(other == veh)
				{
					other = other.next;
					continue;
				}
				
				if( ( approximateDistance = nearby(veh, other, separateDist, separateAngleCos) ) )
				{
					dx = px - other.position.x;
					dy = py - other.position.y;
					dz = pz - other.position.z;
					
					f = 1 / (approximateDistance * approximateDistance);
					
					ax += dx * f;
					ay += dy * f;
					az += dz * f;
				}
				other = other.next;
			}
			
			accumulator.x = ax;
			accumulator.y = ay;
			accumulator.z = az;
			accumulator.setApproximateNormalize();
			
			accumulator.x *= separateStrength;
			accumulator.y *= separateStrength;
			accumulator.z *= separateStrength;
			
			veh.applyGlobalForce(accumulator);
		}
	}
}
