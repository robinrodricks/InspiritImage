package ru.inspirit.fluid.emitter 
{

	/**
	 * @author Eugene Zatepyakin
	 */
	public interface IFluidEmitter 
	{
		function update():void;
		function updateWithForce(fx:Number, fy:Number, fz:Number):void;
		function setRandomPosition():void;
	}
}
