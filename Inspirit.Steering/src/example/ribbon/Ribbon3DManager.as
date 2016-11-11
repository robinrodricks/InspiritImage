package example.ribbon 
{
	import flash.display.Graphics;
	import flash.display.Shape;
	import flash.geom.Matrix3D;
	import flash.geom.PerspectiveProjection;
	import flash.geom.Utils3D;
	import flash.geom.Vector3D;

	/**
	 * @author Eugene Zatepyakin
	 */
	public class Ribbon3DManager 
	{
		protected const projection:PerspectiveProjection = new PerspectiveProjection();
		protected const world:Matrix3D = new Matrix3D();
		protected const transformMatrix:Matrix3D = new Matrix3D();
		
		protected var projectionMatrix:Matrix3D;
		
		public var projected:Vector.<Number>;
		public var uvtData:Vector.<Number>;
		
		public var gfx:Graphics;
		public var viewport:Shape;
		
		public function Ribbon3DManager(viewport:Shape)
		{
			this.viewport = viewport;
			this.gfx = viewport.graphics;
			
			projection.fieldOfView = 60;
			
			projectionMatrix = projection.toMatrix3D();
			
			world.identity(); 
			world.appendTranslation(0, 0, 500);
			world.append(projectionMatrix);
			
			projected = new Vector.<Number>();
			uvtData = new Vector.<Number>();
		}
		
		public function prerender():void
		{
			world.identity(); 
			world.appendRotation(viewport.mouseX/10, Vector3D.Y_AXIS);
			world.appendRotation(-(viewport.mouseY/10), Vector3D.X_AXIS);
			world.appendTranslation(0, 0, 500);
			world.append(projectionMatrix);
			
			gfx.clear();
		}
		
		public function drawRibbon(trans:Matrix3D, color:uint, vertices:Vector.<Number>, indices:Vector.<int>):void
		{
			transformMatrix.identity();
			transformMatrix.append(trans);
			transformMatrix.append(world);
			
			Utils3D.projectVectors( transformMatrix, vertices, projected, uvtData );
			
			gfx.beginFill(color, 0.6);
			gfx.drawTriangles(projected, indices);
			gfx.endFill();
			
			projected.length = uvtData.length = 0;
		}
	}
}
