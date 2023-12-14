function [nS, nD, nN, nB] = categorizeObj(objint_st, cont_st)
    % Separate satellites into satellite, derelict, nonusable/debris, rocket body
    % Input: objint_st: list of object index
    %        cont_st: list of control index
    % Output:   nS: satellite numbers
    %           nD: derelict numbers
    %           nN: debris numbers
    %           nB: rocket body numbers
    % object classes:
    %   P(payload),PMRO(debris),Pfrag(debris),Pdeb(debris),RB(rocket body),
    %   RBMRO(debris),RBfrag(debris),RBdeb(debris),Deb(debris),OtherDeb(debris),
    %   Unkwn(debris),untracked(debris)
    % S: satellite, object index 1 and control index 1
    % D: derelict, object index 1 and control index 0
    % N: nousable/debris, object index 3,4,6,7,8,9..
    % B: rocket body, object index 5
    nS = sum(objint_st==1 & cont_st==1);
    nD = sum(objint_st==1 & cont_st==0);
    nN = sum(objint_st==3 | objint_st==4 | objint_st>=6);
    nB = sum(objint_st==5);
end