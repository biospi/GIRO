classdef (Abstract) Register
    
    properties (Access = protected)
        
        Samples
        
        DeformedSamples
        
        numSamples
        
        RT
        
        sizeRT
        
        sizeDyadicRT
        
        indRT_Start
        
        RT_Adjustment
        
        Ref_Register
        
        Tar_Register
        
    end
    
    methods (Abstract)
        
        OBJ_Register = searching_strategy(Obj_Register)
        
        OBJ_Register = interf_register_solver(OBJ_Register, VAR_Reg)

        % Deformation mechanism:
        OBJ_Register = deforxm_RT(OBJ_Register)
        
    end
    
    methods (Access = public)
                
         % Get an input of a Data object to initiate the Register class:
         function OBJ_Register = Register(OBJ_Data)
            
             OBJ_Register.Samples = OBJ_Data.get_Samples();
             
             OBJ_Register.numSamples = OBJ_Data.get_nuSamples();
             
             [OBJ_Register.RT, OBJ_Register.sizeRT, OBJ_Register.sizeDyadicRT, OBJ_Register.indRT_Start] = OBJ_Data.get_RT_info();
             
             OBJ_Register.RT_Adjustment = zeros(OBJ_Register.numSamples, OBJ_Register.sizeRT);
             
         end
        
         
         % Get the output RT adjustment:
         function [RT, RT_Adjustment] = get_RT_adjustment(OBJ_Register)
        
             RT = OBJ_Register.RT;
             
             RT_Adjustment = OBJ_Register.RT_Adjustment;
             
             
         end
             
         % Get the RT-adjusted samples:
         function DeformedSamples = get_RT_adjusted_samples(OBJ_Register)
             
             DeformedSamples = OBJ_Register.DeformedSamples;
             
         end
             
             
         
    end
    
end