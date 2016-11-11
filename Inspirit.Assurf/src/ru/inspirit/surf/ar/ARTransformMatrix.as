package ru.inspirit.surf.ar 
{

	/**
	 * @author Eugene Zatepyakin
	 */
	public class ARTransformMatrix 
	{
		public var m00:Number;
		public var m01:Number;
		public var m02:Number;
		public var m03:Number;
		public var m10:Number;
		public var m11:Number;
		public var m12:Number;
		public var m13:Number;
		public var m20:Number;
		public var m21:Number;
		public var m22:Number;
		public var m23:Number;
		
		public function clone():ARTransformMatrix
		{
			var m:ARTransformMatrix = new ARTransformMatrix();
			m.m00 = m00;
			m.m01 = m01;
			m.m02 = m02;
			m.m03 = m03;
			m.m10 = m10;
			m.m11 = m11;
			m.m12 = m12;
			m.m13 = m13;
			m.m20 = m20;
			m.m21 = m21;
			m.m22 = m22;
			m.m23 = m23;
			
			return m;
		}

		public function toString():String 
		{
			var str:String = 'transform matrix[3][4]\n';
			str += m00 + '\t' + m01+ '\t' + m02+ '\t' + m03 + '\n';
			str += m10 + '\t' + m11+ '\t' + m12+ '\t' + m13 + '\n';
			str += m20 + '\t' + m21+ '\t' + m22+ '\t' + m23;
			return str;
		}
	}
}
