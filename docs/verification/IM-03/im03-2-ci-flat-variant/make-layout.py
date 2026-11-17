import pathlib
p=pathlib.Path('out/issue102-flat')
for label in ('nested','flat'):
 text=(p/(label+'-probe.cpp')).read_text(encoding='utf-8')
 selected=('std::is_nothrow_constructible_v<FreeStorage,std::in_place_type_t<std::array<double,3>>,const std::array<double,3>&>' if label=='nested' else 'std::is_nothrow_constructible_v<FreeSingleStorage,const FreeSingleStorage&>')
 text+='\nvoid consume(const QueryStorage&){}\nint main(){std::printf("size=%zu selected_nothrow=%d edge_copy=%d edge_move=%d face_copy=%d face_move=%d\\n",sizeof(QueryStorage),int('+selected+'),int(std::is_nothrow_copy_constructible_v<EdgeQueryStorage>),int(std::is_nothrow_move_constructible_v<EdgeQueryStorage>),int(std::is_nothrow_copy_constructible_v<FaceQueryStorage>),int(std::is_nothrow_move_constructible_v<FaceQueryStorage>));}\n'
 (p/(label+'-layout.cpp')).write_bytes(text.encode())
