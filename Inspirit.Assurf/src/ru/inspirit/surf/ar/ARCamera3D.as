package ru.inspirit.surf.ar 
{

	/**
	 * @author Eugene Zatepyakin
	 */
	public class ARCamera3D 
	{		
		// projection matrix
		public var m00:Number = 700.9514702992245;
		public var m01:Number = 0;
		public var m02:Number = 316.5;
		public var m03:Number = 0;
		public var m10:Number = 0;
		public var m11:Number = 726.0941816535367;
		public var m12:Number = 241.5;
		public var m13:Number = 0;
		public var m20:Number = 0;
		public var m21:Number = 0;
		public var m22:Number = 1.0;
		public var m23:Number = 0;
		
		public var screenW:int = 640;
		public var screenH:int = 480;
		
		public function ARCamera3D(screenWidth : int = 640, screenHeight : Number = 480) 
		{
			changeScreenSize(screenWidth, screenHeight);
		}
		
		public final function changeScreenSize(w:int, h:int):void 
		{
			const scale:Number = w / screenW;
			
			this.m00 = this.m00 * scale;
			this.m10 = this.m10 * scale;
			this.m01 = this.m01 * scale;
			this.m11 = this.m11 * scale;
			this.m02 = this.m02 * scale;
			this.m12 = this.m12 * scale;
			this.m03 = this.m03 * scale;
			this.m13 = this.m13 * scale;
			
			screenW = w;
			screenH = h;
		}
		
		public final function identity():void
		{
			m00 = 1.0;
			m01 = 0.0;
			m02 = 0.0;
			m03 = 0.0;
			m10 = 0.0;
			m11 = 1.0;
			m12 = 0.0;
			m13 = 0.0;
			m20 = 0.0;
			m21 = 0.0;
			m22 = 1.0;
			m23 = 0.0;
		}
		
		public final function toString():String 
		{
			var str:String = 'cam matrix[3][4]\t' + screenW + 'x' + screenH + '\n';
			str += m00 + '\t' + m01+ '\t' + m02+ '\t' + m03 + '\n';
			str += m10 + '\t' + m11+ '\t' + m12+ '\t' + m13 + '\n';
			str += m20 + '\t' + m21+ '\t' + m22+ '\t' + m23;
			return str;
		}
	}
}
