classdef UF
    %{
    An implementation of union find data structure.
    It uses weighted quick union by rank with path compression.
    %}
    properties
    %{
    Initialize an empty union find object with N items.
        
    Args:
    	N: Number of items in the union find object.
        
     %}

        id_ % id_ is an array containing N elements, and the ith element value is the ID of the data point contained in the ith point.
        count_
        rank_
    end
    
    methods % Multiple ordinary functions can be placed in the same methods, but if they are static functions, you must place one static function in each method.
        function self = UF(N)
            self.id_ = (1:N)';
            self.count_ = N;
            self.rank_ = zeros(N,1);
        end
        
       
        function [p,self] = find(self, p) 
        % The find function uses a chained list to find the root node, or class number, in the class.

            %Find the set identifier for the item p.
            id = self.id_;
            while (p ~= id(p)) 
               
                id(p) = id(id(p)); 
                p = id(p);
            end
            %  To actually implement this tree shape transformation, 
            % assign the local variable id to the member variable id_.
            self.id_ = id;
        end
        
        
        function num = count(self)
        %Return the number of items.

            num = self.count_;
        end
        
        function iseq = connected(self, p, q)
        % Check if the items p and q are on the same set or not.
        
            iseq = (self.find(p) == self.find(q)); 
        end
        
        
        % When a join is created between two data points, the merging of two classes (or trees) is realized
        % In order to improve the traversal, the main idea of the merge is to add the tree with a shallow level to the tree with a deeper level.
        % i.e., the role of rank is to count the depth of the hierarchy, with larger values indicating deeper trees.
        function self = union(self, p, q)
            % Combine sets containing p and q into a single set.
            
            rank = self.rank_;   
            num = self.count_;  
            [i,self] = self.find(p); % Find the tree that p belongs to and update the tree structure so that all nodes on the path from p to the root node point to the root node
            [j,self] = self.find(q); % Find the tree to which q belongs and update the tree structure so that nodes on the path from q to the root point to the root node          
            id = self.id_;  % Using the updated tree         
            
            if (i == j) % Both data points are in the same category and do not need to be updated
            else
                num = num - 1;
                % When tree i is shallower than tree j, add tree i to tree j, i.e. root node of i points to j
                if (rank(i) < rank(j))
                    id(i) = j; % The root node of tree i points to j
                % When tree j is shallower than tree i, add tree j to tree i, i.e., the root node of j points to i
                elseif(rank(i) > rank(j))
                    id(j) = i; % The root node of the % j tree points to i
                else  % When trees i and j are the same depth, specify that tree j is added to tree i and the root node of j points to i
                    id(j) = i; % The root node of % j points to i.
                    rank(i) = rank(i) + 1; % As a result, the i-tree becomes deeper
                end
            end
            self.id_ = id;
            self.rank_ = rank;
            self.count_ = num;
        end
    end
end
