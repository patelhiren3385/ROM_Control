function obj = get_nodes_coords_connectivity(dim, elements, dof_per_node)
obj.length = dim.length;
obj.width = dim.width;
obj.depth = dim.depth;
%==========================================================================
obj.total_elements = elements.beam;
obj.total_nodes = elements.beam + 1;
obj.total_dofs = obj.total_nodes * dof_per_node;
obj.connectivity = transpose([1:(obj.total_nodes-1);2:obj.total_nodes]);
obj.dof_per_element = dof_per_node * 2;
%==========================================================================
obj.length_sec1 = dim.length_sec1 ; 
obj.elelen_sec1 = dim.length_sec1/elements.sec1 ;
obj.total_elements_sec1 = elements.sec1 ; 
obj.total_nodes_sec1 = elements.sec1 + 1 ;
obj.connectivity_sec1 = transpose([1:(obj.total_nodes_sec1-1);2:obj.total_nodes_sec1]);

obj.length_sec2 = dim.length_sec2 ; 
obj.elelen_sec2 = dim.length_sec2/elements.sec2 ;
obj.total_elements_sec2 = elements.sec2 ; 
obj.total_nodes_sec2 = elements.sec2 + 1 ;
obj.connectivity_sec2 = transpose([1:(obj.total_nodes_sec2-1);2:obj.total_nodes_sec2]);

obj.length_sec3 = dim.length_sec3 ; 
obj.elelen_sec3 = dim.length_sec3/elements.sec3 ;
obj.total_elements_sec3 = elements.sec3 ; 
obj.total_nodes_sec3 = elements.sec3 + 1 ;
obj.connectivity_sec3 = transpose([1:(obj.total_nodes_sec3-1);2:obj.total_nodes_sec3]);

obj.length_sec4 = dim.length_sec4 ; 
obj.total_elements_sec4 = elements.sec4 ; 
obj.elelen_sec4 = dim.length_sec4/elements.sec4 ;
obj.total_nodes_sec4 = elements.sec4 + 1 ;
obj.connectivity_sec4 = transpose([1:(obj.total_nodes_sec4-1);2:obj.total_nodes_sec4]);

obj.length_sec5 = dim.length_sec5 ; 
obj.elelen_sec5 = dim.length_sec5/elements.sec5 ;
obj.total_elements_sec5 = elements.sec5 ; 
obj.total_nodes_sec5 = elements.sec5 + 1 ;
obj.connectivity_sec5 = transpose([1:(obj.total_nodes_sec5-1);2:obj.total_nodes_sec5]);

obj.length_sec6 = dim.length_sec6 ; 
obj.elelen_sec6 = dim.length_sec6/elements.sec6 ;
obj.total_elements_sec6 = elements.sec6 ; 
obj.total_nodes_sec6 = elements.sec6 + 1 ;
obj.connectivity_sec6 = transpose([1:(obj.total_nodes_sec6-1);2:obj.total_nodes_sec6]);

obj.length_sec7 = dim.length_sec7 ; 
obj.elelen_sec7 = dim.length_sec7/elements.sec7 ;
obj.total_elements_sec7 = elements.sec7 ; 
obj.total_nodes_sec7 = elements.sec7 + 1 ;
obj.connectivity_sec7 = transpose([1:(obj.total_nodes_sec7-1);2:obj.total_nodes_sec7]);


%==========================================================================
obj.global_coordinates_sec1 = linspace(0,obj.length_sec1 , obj.total_nodes_sec1); % Fixed Length to Act 1 start
obj.global_coordinates_sec2 = linspace(obj.length_sec1  , obj.length_sec2+obj.length_sec1,obj.total_nodes_sec2); % Act 1 
obj.global_coordinates_sec3 = linspace(obj.length_sec1 + obj.length_sec2  , obj.length_sec3+obj.length_sec2+obj.length_sec1,obj.total_nodes_sec3); % Act 1 to Act 2 start
obj.global_coordinates_sec4 = linspace(obj.length_sec1 + obj.length_sec2 + obj.length_sec3 , obj.length_sec4+obj.length_sec3+obj.length_sec2+obj.length_sec1,obj.total_nodes_sec4); % Act 2 
obj.global_coordinates_sec5 = linspace(obj.length_sec1 + obj.length_sec2 + obj.length_sec3 + obj.length_sec4  , obj.length_sec5 + obj.length_sec4 + obj.length_sec3+obj.length_sec2+obj.length_sec1,obj.total_nodes_sec5); % Act 2 to Act 3 start
obj.global_coordinates_sec6 = linspace(obj.length_sec1 + obj.length_sec2 + obj.length_sec3 + obj.length_sec4 + obj.length_sec5  , obj.length_sec6 + obj.length_sec5 + obj.length_sec4+obj.length_sec3+obj.length_sec2+obj.length_sec1,obj.total_nodes_sec6); %  Act 3 
obj.global_coordinates_sec7 = linspace(obj.length_sec1 + obj.length_sec2 + obj.length_sec3 + obj.length_sec4 + obj.length_sec5 + obj.length_sec6  , obj.length_sec7 + obj.length_sec6 + obj.length_sec5 + obj.length_sec4+obj.length_sec3+obj.length_sec2+obj.length_sec1,obj.total_nodes_sec7); % Act 3 to end
%==========================================================================
obj.element_coord_sec1 = zeros(obj.total_elements_sec1,2);
obj.element_coord_sec2 = zeros(obj.total_elements_sec2,2);
obj.element_coord_sec3 = zeros(obj.total_elements_sec3,2);
obj.element_coord_sec4 = zeros(obj.total_elements_sec4,2);
obj.element_coord_sec5 = zeros(obj.total_elements_sec5,2);
obj.element_coord_sec6 = zeros(obj.total_elements_sec6,2);
obj.element_coord_sec7 = zeros(obj.total_elements_sec7,2);
%==========================================================================
for i=1:obj.total_elements_sec1
    obj.element_coord_sec1(i,:) = obj.global_coordinates_sec1(obj.connectivity_sec1(i,:));
end

for j=1:obj.total_elements_sec2
    obj.element_coord_sec2(j,:) = obj.global_coordinates_sec2(obj.connectivity_sec2(j,:)) ;
end

for i=1:obj.total_elements_sec3
    obj.element_coord_sec3(i,:) = obj.global_coordinates_sec3(obj.connectivity_sec3(i,:));
end

for i=1:obj.total_elements_sec4
    obj.element_coord_sec4(i,:) = obj.global_coordinates_sec4(obj.connectivity_sec4(i,:));
end

for i=1:obj.total_elements_sec5
    obj.element_coord_sec5(i,:) = obj.global_coordinates_sec5(obj.connectivity_sec5(i,:));
end

for i=1:obj.total_elements_sec6
    obj.element_coord_sec6(i,:) = obj.global_coordinates_sec6(obj.connectivity_sec6(i,:));
end

for i=1:obj.total_elements_sec7
    obj.element_coord_sec7(i,:) = obj.global_coordinates_sec7(obj.connectivity_sec7(i,:));
end
%==========================================================================
obj.element_coordinates = [obj.element_coord_sec1 ; obj.element_coord_sec2 ; obj.element_coord_sec3 ; obj.element_coord_sec4 ; obj.element_coord_sec5 ; obj.element_coord_sec6 ; obj.element_coord_sec7 ] ;
%==========================================================================
obj.dofs = node2dof(1:obj.total_nodes, dof_per_node);
if strcmp(dim.support_condition,'c')
    obj.arrested_nodes = 1;
    obj.arrested_dofs = node2dof(obj.arrested_nodes, dof_per_node);
    obj.arrested_dofs = obj.arrested_dofs(:);
    obj.free_dofs = setdiff(1:obj.total_dofs, obj.arrested_dofs);
elseif strcmp(dim.support_condition,'ss')
    obj.arrested_nodes = [1, obj.total_nodes];
    obj.arrested_dofs = node2dof(obj.arrested_nodes, dof_per_node);
    obj.arrested_dofs(end,:) = [];
    obj.arrested_dofs = obj.arrested_dofs(:);
    obj.free_dofs = setdiff(1:obj.total_dofs, obj.arrested_dofs);
end
end