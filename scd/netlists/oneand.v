module oneand(g_input, e_input, o);

	input [0:0] g_input;
	input [0:0] e_input;
	output [0:0] o;

	assign o = g_input & e_input;
	//AND U ( .A(g_input), .B(e_input), .Z(o) );

endmodule