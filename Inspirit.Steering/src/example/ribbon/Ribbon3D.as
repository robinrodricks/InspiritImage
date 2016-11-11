package example.ribbon {
	import ru.inspirit.steering.SteerVector3D;
	import ru.inspirit.steering.Vehicle;
	import ru.inspirit.steering.behavior.BehaviorList;

	import flash.geom.Matrix3D;
	import flash.geom.Vector3D;

	/**
	 * @author Eugene Zatepyakin
	 */
	public final class Ribbon3D extends Vehicle
	{		
		public static var SCALE:Number = 10;
		
		public var x:Number;
		public var y:Number;
		public var z:Number;
		public var length:int;
		
		public var color:uint;
		public var width:Number;
		public var maxLength:uint;
		
		public var vertices:Vector.<Number>;
		public var tvertices:Vector.<Number>;
		public var indices:Vector.<int>;
		
		public var manager:Ribbon3DManager;
		
		public var transform:Matrix3D = new Matrix3D();
		public var target:Vector3D = new Vector3D();
		
		public function Ribbon3D(manager:Ribbon3DManager, position:SteerVector3D = null, behaviorList:BehaviorList = null, color:uint = 0x000000, width:Number = 5, length:uint = 100)
		{
			super(position, behaviorList);
			
			this.manager = manager;
			
			x = 0;
			y = 0;
			z = 0;
			
			this.color = color;
			this.width = width;
			this.maxLength = length;
			
			vertices = new Vector.<Number>();
			tvertices = new Vector.<Number>();
			indices = new Vector.<int>();
			
			vertices.push( x + width, y, z );
			vertices.push( x - width, y, z );
			
			this.length = 2;
			
			var n:int = 4;
			var i:int = maxLength;
			while(--i)
			{
				indices.push(n - 1, n - 2, n - 3);
				indices.push(n - 2, n - 4, n - 3);
				n += 2;
			}
		}

		override public function update():void
		{
			super.update();
			
			target.x = x = position.x * SCALE;
			target.y = y = position.y * SCALE;
			target.z = z = position.z * SCALE;
			
			vertices.push( x + width, y, z );
			vertices.push( x - width, y, z );
			
			length += 2;
			
			while (length > maxLength)
			{
				vertices.shift();
				vertices.shift();
				vertices.shift();
				
				length--;
			}
			
			width += (velocity.length*10 - width) * 0.25;
			
			// You can apply additional transform to each ribbon before render
			transform.identity();
			//transform.appendTranslation(target.x, target.y, target.z);
			//transform.pointAt(target.add(velocity), Vector3D.Z_AXIS, Vector3D.Y_AXIS);
			
			manager.drawRibbon(transform, color, vertices, indices);
		}
	}
}
