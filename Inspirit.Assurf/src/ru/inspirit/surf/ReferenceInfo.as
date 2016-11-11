package ru.inspirit.surf 
{
	import flash.geom.Point;
	import apparat.memory.Memory;

	import ru.inspirit.surf.ar.ARTransformMatrix;

	/**
	 * @author Eugene Zatepyakin
	 */
	public final class ReferenceInfo 
	{
		public var id:int;
		public var width:int;
		public var height:int;
		public var pointsCount:int;
		public var matchedPointsCount:int;
		
		public var transformError:Number;
		
		public const homography:HomographyMatrix = new HomographyMatrix();
		public const transform:ARTransformMatrix = new ARTransformMatrix();
		
		internal var pointsCountPtr:int;
		internal var matchedPointsCountPtr:int;
		internal var transformErrorPtr:int;
		internal var homographyPtr:int;
		internal var transformPtr:int;
		
		internal var pt0:Point;
		internal var pt1:Point;
		internal var pt2:Point;
		internal var pt3:Point;
		
		public function ReferenceInfo(width:int, height:int)
		{
			this.width = width;
			this.height = height;
			
			pt0 = new Point(0, 0);
			pt1 = new Point(width, 0);
			pt2 = new Point(width, height);
			pt3 = new Point(0, height);
		}

		final internal function updateInfo():void
		{
			pointsCount = Memory.readInt(pointsCountPtr);
			matchedPointsCount = Memory.readInt(matchedPointsCountPtr);
			transformError = Memory.readDouble(transformErrorPtr);
		}
		
		final internal function updateHomography():void
		{
			var hp:int = homographyPtr;
			homography.m11 = Memory.readDouble(hp);
			homography.m12 = Memory.readDouble(hp + 8);
			homography.m13 = Memory.readDouble(hp + 16);
			homography.m21 = Memory.readDouble(hp + 24);
			homography.m22 = Memory.readDouble(hp + 32);
			homography.m23 = Memory.readDouble(hp + 40);
			homography.m31 = Memory.readDouble(hp + 48);
			homography.m32 = Memory.readDouble(hp + 56);
			homography.m33 = Memory.readDouble(hp + 64);
		}
		
		final internal function updateTransform():void
		{
			var tp:int = transformPtr;
			transform.m00 = Memory.readDouble(tp);
			transform.m01 = Memory.readDouble(tp + 8);
			transform.m02 = Memory.readDouble(tp + 16);
			transform.m03 = Memory.readDouble(tp + 24);
			transform.m10 = Memory.readDouble(tp + 32);
			transform.m11 = Memory.readDouble(tp + 40);
			transform.m12 = Memory.readDouble(tp + 48);
			transform.m13 = Memory.readDouble(tp + 56);
			transform.m20 = Memory.readDouble(tp + 64);
			transform.m21 = Memory.readDouble(tp + 72);
			transform.m22 = Memory.readDouble(tp + 80);
			transform.m23 = Memory.readDouble(tp + 88);
		}
	}
}
