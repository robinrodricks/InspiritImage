package ru.inspirit.steering 
{

	import flash.geom.Vector3D;
	
	/**
	 * Steer Vector is based on Vector3D
	 * with some additional methods
	 * 
	 * @author Eugene Zatepyakin
	 */
	public final class SteerVector3D extends Vector3D
	{
		public static const ZERO:SteerVector3D = new SteerVector3D();
		
		public static var seed:uint = 1;//int(Math.random() * 0x7FFFFFFE) + 1;
		
		public function SteerVector3D(x:Number = 0, y:Number = 0, z:Number = 0) 
		{
			super(x, y, z, 0);
		}
		
		public function setTo(x:Number = 0, y:Number = 0, z:Number = 0):void
		{
			this.x = x;
			this.y = y;
			this.z = z;
		}
		
		public function setToVector(v:SteerVector3D):void
		{
			this.x = v.x;
			this.y = v.y;
			this.z = v.z;
		}
		
		public function setZero():void
		{
			x = y = z = 0;
		}
		
		public function equalsZero():Boolean
		{
			if (x != 0) return false;
			if (y != 0) return false;
			if (z != 0) return false;
			
			return true;
		}
		
		public function setSum(a:SteerVector3D, b:SteerVector3D):void
		{
			x = a.x + b.x;
			y = a.y + b.y;
			z = a.z + b.z;
		}

		public function scale(factor:Number):SteerVector3D
		{
			return new SteerVector3D(x * factor, y * factor, z * factor);
		}
		
		public function setScale(factor:Number, other:SteerVector3D):void
		{
			x = other.x * factor;
			y = other.y * factor;
			z = other.z * factor;
		}

		public function magnitudeSquared():Number
		{
			return x * x + y * y + z * z;
		}

		public function magnitude():Number
		{
			return Math.sqrt( x * x + y * y + z * z );
		}

		public function truncate(threshold:Number):SteerVector3D
		{
			var lengthSquared:Number = x * x + y * y + z * z;
			
			if(lengthSquared > threshold * threshold)
			{
				var factor:Number = threshold / Math.sqrt(lengthSquared);
				return new SteerVector3D(x * factor, y * factor, z * factor);
			} 
			else
			{
				return this;
			}
		}
		
		public function setApproximateTruncate(threshold:Number):void
		{
			var length:Number = approximateLength();
			if (length > threshold) {
				length = threshold / length;
				x *= length;
				y *= length;
				z *= length;
			}
		}
		
		public function approximateTruncate(threshold:Number):SteerVector3D
		{
			var length:Number = approximateLength();
			if(length > threshold)
			{
				var factor:Number = threshold / length;
				return new SteerVector3D(x * factor, y * factor, z * factor);
			}
			else
			{
				return this;
			}
		}

		public function dot(that:SteerVector3D):Number
		{
			return x * that.x + y * that.y + z * that.z;
		}

		public function cross(that:SteerVector3D):SteerVector3D
		{
			return new SteerVector3D(y * that.z - z * that.y, z * that.x - x * that.z, x * that.y - y * that.x);
		}
		
		public function setDiff(a:SteerVector3D, b:SteerVector3D):void
		{
			x = a.x - b.x;
			y = a.y - b.y;
			z = a.z - b.z;
		}
		
		public function setCross(a:SteerVector3D, b:SteerVector3D):void
		{
			x = a.y * b.z - a.z * b.y;
			y = a.z * b.x - a.x * b.z;
			z = a.x * b.y - a.y * b.x;
		}
		
		public function setUnitRandom():void
		{
			do
			{
				x = random() * 2 - 1;
				y = random() * 2 - 1;
				z = random() * 2 - 1;
			} while((x * x + y * y + z * z) > 1);
		}

		public static function unitRandom():SteerVector3D
		{
			var x:Number = 0, y:Number = 0, z:Number = 0;
			do
			{
				x = random() * 2 - 1;
				y = random() * 2 - 1;
				z = random() * 2 - 1;
			} while((x * x + y * y + z * z) > 1);
			
			return new SteerVector3D(x, y, z);
		}

		public function approximateLength():Number
		{
			var a:Number = x;
			if (a < 0) a = -a;
			
			var b:Number = y;
			if (b < 0) b = -b;
			
			var c:Number = z;
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
			return a * 0.9375 + (b + c) * 0.375;
		}
		
		public function distance(v:SteerVector3D):Number
		{
			var xd:Number = x - v.x;
			var yd:Number = y - v.y;
			var zd:Number = z - v.z;
			
			return Math.sqrt( xd * xd + yd * yd + zd * zd );
		}
		
		public function approximateDistance(v:SteerVector3D):Number
		{
			var xd:Number = x - v.x;
			var yd:Number = y - v.y;
			var zd:Number = z - v.z;
			
			var a:Number = xd;
			if (a < 0) a = -a;
			
			var b:Number = yd;
			if (b < 0) b = -b;
			
			var c:Number = zd;
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
			return a * 0.9375 + (b + c) * 0.375;
		}
		
		public function setApproximateNormalize():void
		{
			var m:Number = approximateLength();
			if(m != 0) 
			{
				m = 1 / m;
				x *= m;
				y *= m;
				z *= m;
			}
		}

		public function approximateNormalize():SteerVector3D
		{
			var m:Number = approximateLength();
			if(m == 0)
			{
				return this;
			}
			else
			{		
				m = 1 / m;
				return new SteerVector3D(x * m, y * m, z * m);
			}
		}
		
		public function setInterp(blend:Number, v0:SteerVector3D, v1:SteerVector3D):void
		{
			x = v0.x + blend * (v1.x - v0.x);
			y = v0.y + blend * (v1.y - v0.y);
			z = v0.z + blend * (v1.z - v0.z);
		}
		
		 public function setNormalize():void
		{
			var m:Number = Math.sqrt( x * x + y * y + z * z );
			if(m != 0) 
			{
				m = 1 / m;
				x *= m;
				y *= m;
				z *= m;
			}
		}
		
		public function parallelComponent (unitBasis:SteerVector3D):SteerVector3D
        {
            const projection:Number = dot(unitBasis);
            return unitBasis.scale(projection);
        }

        // return component of vector perpendicular to a unit basis vector
        // (IMPORTANT NOTE: assumes "basis" has unit magnitude (length==1))

        /*public function perpendicularComponent (unitBasis:SteerVector3D):SteerVector3D
        {
            return subtract(parallelComponent(unitBasis));
        }*/
        
        public static function random():Number
        {
        	return (seed = (seed * 16807) % 2147483647) / 2147483647;
        }
	}
}
